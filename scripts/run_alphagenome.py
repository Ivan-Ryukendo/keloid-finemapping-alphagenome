"""Score keloid variants with AlphaGenome 2 (v2 API) in skin cell-type conditions.

Input : data/raw/keloid_gwas_variants.csv
Output: results/alphagenome/variant_scores_long.csv   (one row per variant x condition)
        results/alphagenome/variant_skipped.csv       (variants not scored, with reason)
Effect metric: DIFF_LOG2_SUM for chromatin/TF/CAGE tracks; GeneMaskLFC (RNA-seq) for genes.
"""
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client as c, variant_scorers as vs

ROOT = Path(__file__).resolve().parent.parent
load_dotenv(ROOT / ".env")
IN = ROOT / "data" / "raw" / "keloid_gwas_variants.csv"
OUTDIR = ROOT / "results" / "alphagenome"
OUTDIR.mkdir(parents=True, exist_ok=True)
LONG = OUTDIR / "variant_scores_long.csv"

CONDITIONS = {  # verified present in supported_conditions()
    "skin_fibroblast": "UBERON:0002097_CL:0002620",
    "dermal_fibroblast": "UBERON:0002097_CL:0002551",
    "myofibroblast": "UBERON:0002097_CL:0000186",
    "keratinocyte": "UBERON:0002097_CL:0000312",
    "macrophage": "UBERON:0002097_CL:0000235",
    "endothelial": "UBERON:0002097_CL:0000115",
}
HALF = 262144  # interval = 524288 bp = SEQUENCE_LENGTH_500KB
BASES = set("ACGT")


def ref_bases(variants):
    """GRCh38 reference base per (chrom,pos) via Ensembl batch endpoint."""
    out = {}
    keys = sorted(set(variants))
    for i in range(0, len(keys), 40):
        chunk = keys[i:i + 40]
        regs = [f"{ch}:{p}..{p}:1" for ch, p in chunk]
        r = requests.post("https://rest.ensembl.org/sequence/region/human",
                          json={"regions": regs}, headers={"Accept": "application/json"}, timeout=60)
        r.raise_for_status()
        for k, item in zip(chunk, r.json()):
            out[k] = item["seq"].upper()
    return out


def prepare():
    d = pd.read_csv(IN)
    # position-only IDs like "chr15:56140995" have no rsID annotation
    m = d.rs_id.str.extract(r"^chr(\w+):(\d+)$")
    d["chrom"] = d.chrom.astype("string").fillna(m[0])
    d["pos_grch38"] = pd.to_numeric(d.pos_grch38, errors="coerce").fillna(pd.to_numeric(m[1]))
    d = d.drop_duplicates("rs_id").reset_index(drop=True)
    refs = ref_bases(zip(d.chrom.astype(str), d.pos_grch38.astype(int)))
    rows, skipped = [], []
    for r in d.itertuples():
        ch, pos = str(r.chrom), int(r.pos_grch38)
        ref = refs[(ch, pos)]
        eff = r.effect_allele if isinstance(r.effect_allele, str) else None
        known = set(re.findall(r"\b[ACGT]\b", str(r.alleles).split("(")[0])) if isinstance(r.alleles, str) else set()
        known = {a for a in known if a in BASES}
        if ref not in BASES:
            skipped.append((r.rs_id, f"ref base {ref!r}")); continue
        others = sorted(known - {ref})
        if eff in BASES and eff != ref:
            alts, role = [eff], "effect_allele_is_alt"
        elif len(others) == 1:
            alts, role = others, ("effect_allele_is_ref" if eff == ref else "effect_allele_unknown")
        elif others:  # multiallelic and the effect allele is the reference: score every alt
            alts, role = others, "multiallelic_all_alts"
        else:
            skipped.append((r.rs_id, f"no alt allele known (ref={ref}, effect={eff})")); continue
        for alt in alts:
            rows.append(dict(rs_id=r.rs_id, chrom=f"chr{ch}", pos=pos, ref=ref, alt=alt, allele_role=role,
                             gene=r.mapped_genes, consequence=r.most_severe_consequence,
                             in_catalog=bool(r.in_gwas_catalog_keloid)))
    pd.DataFrame(skipped, columns=["rs_id", "reason"]).to_csv(OUTDIR / "variant_skipped.csv", index=False)
    return pd.DataFrame(rows)


def prepare_finemap():
    """LD-proxy SNVs (finemap_candidates.csv); ref/alt already GRCh38, ref checked against Ensembl."""
    d = pd.read_csv(ROOT / "data" / "processed" / "finemap_candidates.csv")
    d = d.drop_duplicates(["proxy_rsid", "alt"]).reset_index(drop=True)
    refs = ref_bases(zip(d.chrom.astype(str), d.pos_grch38.astype(int)))
    bad = [r.proxy_rsid for r in d.itertuples() if refs[(str(r.chrom), int(r.pos_grch38))] != r.ref]
    assert not bad, f"ref mismatch: {bad[:5]}"
    return pd.DataFrame(dict(rs_id=d.proxy_rsid, chrom="chr" + d.chrom.astype(str), pos=d.pos_grch38.astype(int),
                             ref=d.ref, alt=d.alt, allele_role="ld_proxy", gene=d.lead_rsid,
                             consequence=d.consequence, in_catalog=False))


def summarize(results):
    """Reduce the scorer AnnData list to one dict of effect sizes."""
    out = {}
    for a in results:
        s = str(a.uns["variant_scorer"])
        if a.shape[0] == 0 or a.shape[1] == 0:
            continue
        x = np.asarray(a.X, dtype=float)
        if "DIFF_LOG2_SUM" in s:
            name = re.search(r"requested_output=(\w+)", s).group(1)
            out[f"{name}_log2fc_absmax"] = float(np.nanmax(np.abs(x)))
        elif s.startswith("GeneMaskLFCScorer(requested_output=RNA_SEQ"):
            i, j = np.unravel_index(np.nanargmax(np.abs(x)), x.shape)
            out["RNA_lfc_absmax"] = float(x[i, j])
            out["RNA_top_gene"] = a.obs["gene_name"].iloc[i]
        elif s.startswith("GeneMaskLFCScorer(requested_output=SC_GENE"):
            out["scRNA_lfc_absmax"] = float(np.nanmax(np.abs(x)))
    return out


def main():
    cl = c.create(os.environ["ALPHAGENOME_API_KEY"])
    avail = set(cl.supported_conditions().condition)
    assert all(v in avail for v in CONDITIONS.values()), "a configured condition is not supported"
    scorers = vs.get_recommended_scorers(c.Organism.HOMO_SAPIENS)
    global LONG
    if "--finemap" in sys.argv:
        LONG = OUTDIR / "finemap_scores_long.csv"
        var = prepare_finemap()
    else:
        var = prepare()
    print(f"{len(var)} variants to score x {len(CONDITIONS)} conditions; skipped list in variant_skipped.csv")
    done = pd.read_csv(LONG) if LONG.exists() else pd.DataFrame()
    seen = set(zip(done.get("rs_id", []), done.get("alt", []), done.get("condition_name", [])))
    rows = done.to_dict("records")
    for v in var.itertuples():
        for cname, cid in CONDITIONS.items():
            if (v.rs_id, v.alt, cname) in seen:
                continue
            itv = genome.Interval(v.chrom, max(0, v.pos - HALF), v.pos + HALF)
            var_obj = genome.Variant(v.chrom, v.pos, v.ref, v.alt)
            for attempt in range(4):
                try:
                    res = cl.score_variant(itv, var_obj, condition=cid, variant_scorers=scorers)
                    break
                except Exception as e:  # rate limits / transient errors
                    print(f"  retry {v.rs_id} {cname}: {type(e).__name__} {str(e)[:80]}", file=sys.stderr)
                    time.sleep(5 * (attempt + 1))
            else:
                print(f"  FAILED {v.rs_id} {cname}", file=sys.stderr); continue
            rows.append({**v._asdict(), "condition_name": cname, "condition_id": cid, **summarize(res)})
            pd.DataFrame(rows).to_csv(LONG, index=False)  # resumable
        print(v.rs_id, v.alt, "done", flush=True)
    out = pd.DataFrame(rows)
    assert len(out) == len(var) * len(CONDITIONS), "not all variant x condition pairs scored"
    print("rows:", len(out))


if __name__ == "__main__":
    main()
