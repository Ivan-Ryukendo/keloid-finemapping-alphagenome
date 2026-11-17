"""Concordance of AlphaGenome 2 keloid-variant scores with PRIOR (independent) evidence.

No gold standard exists for keloid causal variants, so this reports CONCORDANCE, not accuracy.
Prior-evidence labels use only GWAS p-values / replication / literature seeds / GTEx eQTL / EuropePMC,
never AlphaGenome output.

Usage (from workspace root, with venv python):
  python scripts/benchmark_vs_prior.py --fetch avi|gtex|ucsc|vep|epmc   # network stages, cached in cache/benchmark/
  python scripts/benchmark_vs_prior.py                                  # analysis from cache (writes CSV + report)

Network stages call: AlphaGenome Atlas API (documented SDK, key from .env, never printed); the project
skill wrappers for GTEx, UCSC (phyloP/phastCons, ENCODE cCRE track), Ensembl VEP, Europe PMC via `uv run`.
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SK = ROOT.parent / ".claude" / "skills"
CACHE = ROOT / "cache" / "benchmark"
CACHE.mkdir(parents=True, exist_ok=True)
GW = 5e-8
RNG = np.random.default_rng(20261116)
B = 5000  # bootstrap / permutation replicates
GTEX_TISSUES = ["Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Cells_Cultured_fibroblasts"]
METRICS = ["ATAC_log2fc_absmax", "CAGE_log2fc_absmax", "CHIP_HISTONE_log2fc_absmax", "CHIP_TF_log2fc_absmax",
           "DNASE_log2fc_absmax", "SC_ATAC_log2fc_absmax", "RNA_lfc_absmax", "scRNA_lfc_absmax"]
PROVIDED7 = METRICS[:7]  # ad hoc ranking used these 7 (RNA as |lfc|)
CHROMATIN = ["ATAC_log2fc_absmax", "CHIP_HISTONE_log2fc_absmax", "CHIP_TF_log2fc_absmax", "DNASE_log2fc_absmax"]
EXPR = ["CAGE_log2fc_absmax", "RNA_lfc_absmax", "scRNA_lfc_absmax"]


def variants():
    l = pd.read_csv(ROOT / "results/alphagenome/variant_scores_long.csv")
    return l, l.drop_duplicates(["rs_id", "alt"])[["rs_id", "chrom", "pos", "ref", "alt"]].reset_index(drop=True)


def uv(script, *args):
    env = dict(os.environ, UV_CACHE_DIR=str(ROOT / "cache" / "uv"))
    return subprocess.run(["uv", "run", str(SK / script), *map(str, args)], env=env, capture_output=True, text=True)


# ------------------------------------------------------------------ fetch stages
def fetch_avi():
    out = CACHE / "avi.csv"
    from dotenv import load_dotenv
    load_dotenv(ROOT / ".env")  # key stays in env, never printed
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    client = atlas.create(os.environ["ALPHAGENOME_API_KEY"])
    _, v = variants()
    rows = []
    for r in v.itertuples():
        var = genome.Variant.from_str(f"{r.chrom}:{r.pos}:{r.ref}>{r.alt}")
        try:
            a = client.query_variant(var, requested_scorers=["AVI_SCORE"])["AVI_SCORE"]
            raw = float(np.ravel(a.X)[0])
            cdf = float(np.ravel(a.layers["quantiles"])[0])
            tail = max(1e-7, 1 - cdf)
            rows.append((r.rs_id, r.alt, raw, cdf, -10 * math.log10(tail), ""))
        except Exception as e:  # keep the failure visible, do not hide it
            rows.append((r.rs_id, r.alt, np.nan, np.nan, np.nan, type(e).__name__ + ":" + str(e)[:150]))
    pd.DataFrame(rows, columns=["rs_id", "alt", "avi_raw", "avi_cdf", "avi_phred", "error"]).to_csv(out, index=False)


def mapped_genes():
    g = pd.read_csv(ROOT / "data/raw/keloid_gwas_variants.csv")
    m = {}
    for rs, d in g.groupby("rs_id", sort=False):
        m[rs] = sorted({x for v in d.mapped_genes.dropna() for x in str(v).split(";")})
    return m


def fetch_gtex():
    # NOTE: the skill's get-eqtls-in-region (singleTissueEqtlByLocation) returned [] for 100 kb windows that
    # contain known NEDD4 skin eQTLs, so it is unreliable; use per-gene eQTLs (get-gene-eqtls) and match by position.
    for sym in sorted({x for v in mapped_genes().values() for x in v}):
        f = CACHE / f"gtexgene_{sym}.json"
        if not f.exists():
            p = uv("gtex_database/scripts/gtex_cli.py", "resolve-gencode-id", sym, "--output", f)
            if p.returncode:
                f.write_text(json.dumps({"unresolved": True}))
        gid = json.load(open(f)).get("gencode_id")
        e = CACHE / f"gtexeqtl_{sym}.json"
        if gid and not e.exists():
            p = uv("gtex_database/scripts/gtex_cli.py", "get-gene-eqtls", gid, "--tissues", ",".join(GTEX_TISSUES), "--output", e)
            if p.returncode:
                print("GTEX FAIL", sym, p.stderr[-200:])


def fetch_ucsc():
    _, v = variants()
    v = v.drop_duplicates("rs_id")
    f = CACHE / "phylop.json"
    if not f.exists():
        p = uv("ucsc_conservation_and_tfbs/scripts/get_conservation.py", "--coordinates",
               *[f"{r.chrom}:{r.pos}" for r in v.itertuples()], "--output", f)
        if p.returncode:
            print("UCSC cons FAIL", p.stderr[-300:])
    # ENCODE cCRE registry via the UCSC track (SCREEN GraphQL coordinate search returned empty for a
    # known-positive HBB test region, so it is not used)
    f = CACHE / "ccre.json"
    if not f.exists():
        p = uv("ucsc_conservation_and_tfbs/scripts/get_tfbs.py", "--coordinates",
               *[f"{r.chrom}:{r.pos - 5}-{r.pos + 5}" for r in v.itertuples()],
               "--tracks", "encodeCcreCombined", "--output", f)
        if p.returncode:
            print("UCSC cCRE FAIL", p.stderr[-300:])


def fetch_vep():
    g = pd.read_csv(ROOT / "data/raw/keloid_gwas_variants.csv").drop_duplicates("rs_id")
    for rs in g[g.is_coding_like & g.rs_id.str.startswith("rs")].rs_id:
        f = CACHE / f"vep_{rs}.json"
        if not f.exists():
            p = uv("ensembl_database/scripts/ensembl_api.py", "vep", rs, "--species", "human", "--output", f)
            if p.returncode:
                print("VEP FAIL", rs, p.stderr[-200:])


def fetch_epmc():
    g = pd.read_csv(ROOT / "data/raw/keloid_gwas_variants.csv").drop_duplicates("rs_id")
    qs = {rs: f'"{rs}" AND keloid' for rs in g.rs_id if rs.startswith("rs")}
    qs["_zhu2013"] = "AUTH:Zhu AND keloid AND genome-wide AND PUB_YEAR:2013"
    qs["_nakashima2010"] = "EXT_ID:20711176 AND SRC:MED"
    for k, q in qs.items():
        f = CACHE / f"epmc_{k}.json"
        if not f.exists():
            p = uv("literature_search_europepmc/scripts/europepmc_api.py", "search", q, "--max_results", 25,
                   "--result_type", "lite", "--output", f)
            if p.returncode:
                print("EPMC FAIL", k, p.stderr[-200:])


# ------------------------------------------------------------------ prior-evidence labels (no AlphaGenome inputs)
def prior_table():
    g = pd.read_csv(ROOT / "data/raw/keloid_gwas_variants.csv")
    rows = []
    for rs, d in g.groupby("rs_id", sort=False):
        d = d.dropna(subset=["p_value"])
        gws = d[d.p_value < GW]
        seed = g[(g.rs_id == rs)].seed_source.dropna()
        rows.append(dict(rs_id=rs, in_catalog=len(d) > 0,
                         min_p=d.p_value.min() if len(d) else np.nan,
                         n_pmids_gws=gws.pubmed_id.nunique(),
                         seed=len(seed) > 0,
                         seed_classic=bool(seed.str.contains("Nakashima|Zhu", regex=True).any()) if len(seed) else False,
                         coding_like=bool(g[g.rs_id == rs].is_coding_like.iloc[0])))
    t = pd.DataFrame(rows).set_index("rs_id")
    t["neglogp"] = -np.log10(t.min_p)
    t["lab_gws"] = t.min_p < GW                      # genome-wide significant in >=1 study
    t["lab_repl"] = t.n_pmids_gws >= 2               # GW-significant in >=2 distinct publications
    t["lab_repl3"] = t.n_pmids_gws >= 3
    # GTEx skin / fibroblast eQTL at the exact SNP, for the mapped gene(s) only (genes without a GTEx
    # record, or variants with no mapped gene, are "not assessed" = NaN)
    mg = mapped_genes()
    pos = variants()[1].drop_duplicates("rs_id").set_index("rs_id").pos
    ge, gn, ga = [], [], []
    for rs in t.index:
        genes, hits, assessed = mg.get(rs, []), set(), False
        for sym in genes:
            f = CACHE / f"gtexeqtl_{sym}.json"
            if f.exists():
                assessed = True
                for e in json.load(open(f)):
                    if int(e["pos"]) == int(pos[rs]):
                        hits.add(f"{sym}:{e['tissueSiteDetailId']}")
        ge.append(len(hits) > 0 if assessed else np.nan)
        gn.append(";".join(sorted(hits)))
        ga.append(assessed)
    t["lab_gtex"], t["gtex_hits"], t["gtex_assessed"] = ge, gn, ga
    # Europe PMC open-access co-mention of rsID and keloid (OA-only => undercount)
    t["epmc_hits"] = [json.load(open(CACHE / f"epmc_{rs}.json")).get("hitCount", np.nan)
                      if (CACHE / f"epmc_{rs}.json").exists() else np.nan for rs in t.index]
    t["lab_lit"] = t.epmc_hits >= 3  # >=3 OA papers co-mention rsID and keloid (1 = just a GWAS table/supplement)
    return t


# ------------------------------------------------------------------ AlphaGenome recipes
def pr(x):  # percentile rank in (0,1]
    return stats.rankdata(x) / len(x)


def recipes(long, avi):
    long = long.copy()
    ids = list(dict.fromkeys(long.rs_id))
    # max over alts and conditions, per metric
    allc = long.groupby("rs_id")[METRICS].max().loc[ids].fillna(0)  # unscored gene-mask RNA = no predicted effect
    R = {}
    ra = lambda cols, df: pd.DataFrame({c: pr(df[c].values) for c in cols}, index=df.index)
    prov = pd.read_csv(ROOT / "results/alphagenome/variant_ranking.csv")
    R["provided_rank_score (max over alts)"] = prov.groupby("rs_id").rank_score.max().reindex(ids)
    R["provided_rank_score (mean over alts)"] = prov.groupby("rs_id").rank_score.mean().reindex(ids)
    R["rank-average, 7 provided metrics (recomputed)"] = ra(PROVIDED7, allc).mean(axis=1)
    R["rank-average, 8 metrics"] = ra(METRICS, allc).mean(axis=1)
    R["max-of-ranks, 7 metrics"] = ra(PROVIDED7, allc).max(axis=1)
    R["rank-average, chromatin tracks only"] = ra(CHROMATIN, allc).mean(axis=1)
    R["rank-average, expression tracks only"] = ra(EXPR, allc).mean(axis=1)
    for m in METRICS:
        R["single track: " + m.replace("_absmax", "")] = allc[m]
    for cond, d in long.groupby("condition_name"):
        c = d.groupby("rs_id")[METRICS].max().loc[ids].fillna(0)
        R[f"rank-average 7 metrics, {cond} only"] = ra(PROVIDED7, c).mean(axis=1)
    if avi is not None:
        a = avi.groupby("rs_id")
        R["AVI phred (Atlas, max over alts)"] = a.avi_phred.max().reindex(ids)
        R["AVI raw (Atlas, max over alts)"] = a.avi_raw.max().reindex(ids)
    return pd.DataFrame(R)


# ------------------------------------------------------------------ statistics
def auc(score, lab):
    s, l = np.asarray(score, float), np.asarray(lab, bool)
    if l.sum() == 0 or (~l).sum() == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(s[l], s[~l], alternative="two-sided")
    return u / (l.sum() * (~l).sum()), p


def fast_auc(s, l):
    n1 = l.sum()
    if n1 == 0 or n1 == len(l):
        return np.nan
    r = stats.rankdata(s)
    return (r[l].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(l) - n1))


def fast_rho(x, y):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return np.corrcoef(rx, ry)[0, 1]


def boot(fn, n, B=B):
    """percentile bootstrap CI of fn(idx) over resampled variant indices"""
    out = []
    for _ in range(B):
        v = fn(RNG.integers(0, n, n))
        if v == v:
            out.append(v)
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if len(out) > B * 0.5 else (np.nan, np.nan)


def spearman_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    rho, p = stats.spearmanr(x, y)
    lo, hi = boot(lambda i: fast_rho(x[i], y[i]), len(x))
    return rho, lo, hi, p, len(x)


def auc_ci(score, lab):
    s, l = np.asarray(score, float), np.asarray(lab, bool)
    ok = ~np.isnan(s)
    s, l = s[ok], l[ok]
    a, p = auc(s, l)
    lo, hi = boot(lambda i: fast_auc(s[i], l[i]), len(s))
    return a, lo, hi, p, int(l.sum()), int((~l).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", choices=["avi", "gtex", "ucsc", "vep", "epmc"])
    a = ap.parse_args()
    if a.fetch:
        {"avi": fetch_avi, "gtex": fetch_gtex, "ucsc": fetch_ucsc, "vep": fetch_vep, "epmc": fetch_epmc}[a.fetch]()
        return
    analyse()


def analyse():
    long, vtab = variants()
    avi = pd.read_csv(CACHE / "avi.csv") if (CACHE / "avi.csv").exists() else None
    t = prior_table()
    R = recipes(long, avi)
    ids = list(R.index)
    t = t.loc[ids]
    # ---- older predictors
    ph = json.load(open(CACHE / "phylop.json")) if (CACHE / "phylop.json").exists() else {}
    cc = json.load(open(CACHE / "ccre.json")) if (CACHE / "ccre.json").exists() else {}
    pos = vtab.drop_duplicates("rs_id").set_index("rs_id")
    phy, pcons, ccre_any, ccre_reg, ccre_cls = [], [], [], [], []
    for rs in ids:
        ch, p = pos.loc[rs, "chrom"], int(pos.loc[rs, "pos"])
        e = ph.get(f"{ch}:{p}", {}).get("tracks", {})
        phy.append(e.get("phyloP100way", [{}])[0].get("value", np.nan) if e.get("phyloP100way") else np.nan)
        pcons.append(e.get("phastCons100way", [{}])[0].get("value", np.nan) if e.get("phastCons100way") else np.nan)
        items = []
        for k, v in cc.items():
            if k.startswith(f"{ch}:"):
                m = re.match(r".*:(\d+)-(\d+)", k)
                if m and int(m.group(1)) == p - 5:
                    items = v["tracks"].get("encodeCcreCombined", [])
                    items = items if isinstance(items, list) else []
        hit = [i for i in items if i["chromStart"] <= p - 1 < i["chromEnd"]]
        ccre_any.append(bool(hit))
        ccre_cls.append(hit[0]["encodeLabel"] if hit else "")
        ccre_reg.append(bool(hit) and hit[0]["encodeLabel"] in ("PLS", "pELS", "dELS", "DNase-H3K4me3"))
    t["phyloP100"], t["phastCons100"] = phy, pcons
    t["ccre_overlap"], t["ccre_class"], t["ccre_regulatory"] = ccre_any, ccre_cls, ccre_reg
    older_ok = {"phyloP100": not np.all(np.isnan(phy)), "ccre": bool(cc)}
    # coding predictors
    cod = []
    for rs in t.index[t.coding_like]:
        f = CACHE / f"vep_{rs}.json"
        if not f.exists():
            cod.append(dict(rs_id=rs, note="VEP not fetched"))
            continue
        d = json.load(open(f))[0]
        tc = [x for x in d.get("transcript_consequences", []) if "missense_variant" in x.get("consequence_terms", [])]
        sift = [x["sift_score"] for x in tc if "sift_score" in x]
        poly = [x["polyphen_score"] for x in tc if "polyphen_score" in x]
        am = [x["am_pathogenicity"] for x in tc if "am_pathogenicity" in x] + \
             [x["alphamissense"]["am_pathogenicity"] for x in tc if "alphamissense" in x]
        cod.append(dict(rs_id=rs, sift_min=min(sift) if sift else np.nan, sift_call=next((x.get("sift_prediction") for x in tc if "sift_prediction" in x), ""),
                        polyphen_max=max(poly) if poly else np.nan, polyphen_call=next((x.get("polyphen_prediction") for x in tc if "polyphen_prediction" in x), ""),
                        alphamissense_max=max(am) if am else np.nan))
    cod = pd.DataFrame(cod).set_index("rs_id") if cod else pd.DataFrame()

    # ---- main concordance table
    rows = []

    def add(section, name, outcome, est, lo, hi, p, n, npos=np.nan, nneg=np.nan, note=""):
        rows.append(dict(section=section, predictor_or_recipe=name, outcome=outcome, estimate=est, ci_lo=lo, ci_hi=hi,
                         p_value=p, n=n, n_pos=npos, n_neg=nneg, note=note))

    incat = t.in_catalog.values
    labs = {
        "GW-significant in >=1 study (vs in-catalog p>=5e-8)": (t.lab_gws, incat),
        "replicated: GW-sig in >=2 publications (vs other in-catalog)": (t.lab_repl, incat),
        "replicated in >=3 publications (vs other in-catalog)": (t.lab_repl3, incat),
        "literature seed (Nakashima/Zhu/other) (vs non-seed, all 50)": (t.seed, np.ones(len(t), bool)),
        "classic seed Nakashima 2010/Zhu 2013 (vs all other 50)": (t.seed_classic, np.ones(len(t), bool)),
        "GTEx skin/fibroblast eQTL at SNP for mapped gene (vs none, assessed only)": (t.lab_gtex.fillna(False).astype(bool), t.gtex_assessed.values),
        "EuropePMC >=3 OA rsID+keloid papers (vs fewer, all 50)": (t.lab_lit, np.ones(len(t), bool)),
    }
    preds = R.copy()
    preds["OLDER: phyloP100way"] = t.phyloP100.values
    preds["OLDER: phastCons100way"] = t.phastCons100.values
    preds["OLDER: ENCODE cCRE regulatory-class overlap"] = t.ccre_regulatory.astype(float).values
    preds["OLDER: ENCODE cCRE any overlap"] = t.ccre_overlap.astype(float).values
    for name in preds.columns:
        x = preds[name].values
        if np.all(np.isnan(x)) or (name.startswith("OLDER: ENCODE") and not cc):
            continue  # stage not fetched
        sec = "older" if name.startswith("OLDER") else "alphagenome"
        rho, lo, hi, p, n = spearman_ci(x[incat], t.neglogp.values[incat])
        add(sec, name, "Spearman vs -log10(min GWAS p), in-catalog", rho, lo, hi, p, n)
        for ln, (lab, mask) in labs.items():
            lab = lab.values.astype(bool)
            if lab[mask].sum() < 2 or (~lab[mask]).sum() < 2:
                continue
            A, lo, hi, p, npos, nneg = auc_ci(x[mask], lab[mask])
            add(sec, name, "AUC: " + ln, A, lo, hi, p, npos + nneg, npos, nneg)

    # ---- family-wise permutation test across the AlphaGenome recipes (guards against cherry-picking)
    ag = [c for c in R.columns]
    obs_rho = {c: stats.spearmanr(R[c].values[incat], t.neglogp.values[incat])[0] for c in ag}
    best = max(ag, key=lambda c: obs_rho[c])
    worst = min(ag, key=lambda c: obs_rho[c])
    nlp = t.neglogp.values[incat]
    Xr = np.apply_along_axis(stats.rankdata, 0, R[ag].values[incat])
    Xr = (Xr - Xr.mean(0)) / Xr.std(0)
    mx = np.empty(B)
    for b in range(B):
        pr_ = stats.rankdata(RNG.permutation(nlp))
        mx[b] = (Xr * ((pr_ - pr_.mean()) / pr_.std())[:, None]).mean(0).max()
    fw_p = (np.sum(mx >= obs_rho[best]) + 1) / (B + 1)
    add("familywise", f"best recipe: {best}", "max Spearman over all recipes vs permuted -log10p", obs_rho[best],
        np.nan, np.nan, fw_p, int(incat.sum()), note="one-sided family-wise permutation p; null = max rho of all recipes under shuffled p")
    rho_vals = np.array(list(obs_rho.values()))
    add("familywise", "all AG recipes", "median / min / max Spearman across recipes", float(np.median(rho_vals)),
        float(rho_vals.min()), float(rho_vals.max()), np.nan, len(ag), note=f"n recipes={len(ag)}; best={best}; worst={worst}")

    # family-wise permutation for the replicated-label AUC (best of all recipes vs shuffled labels)
    lab = t.lab_repl.values[incat].astype(bool)
    aucs = {c: fast_auc(R[c].values[incat], lab) for c in ag}
    bestA = max(ag, key=lambda c: aucs[c])
    ranks = np.apply_along_axis(stats.rankdata, 0, R[ag].values[incat])
    n1 = lab.sum()
    mxa = np.empty(B)
    for b in range(B):
        pl = RNG.permutation(lab)
        mxa[b] = ((ranks[pl].sum(0) - n1 * (n1 + 1) / 2) / (n1 * (len(lab) - n1))).max()
    fwa = (np.sum(mxa >= aucs[bestA]) + 1) / (B + 1)
    av = np.array(list(aucs.values()))
    add("familywise", f"best recipe: {bestA}", "max AUC (replicated label) over all recipes vs shuffled labels", aucs[bestA], np.nan, np.nan, fwa,
        int(incat.sum()), int(n1), int(len(lab) - n1), note="family-wise permutation p")
    add("familywise", "all AG recipes", "median / min / max AUC (replicated label) across recipes", float(np.median(av)), float(av.min()),
        float(av.max()), np.nan, len(ag))
    # leave-one-variant-out stability of the replicated-label AUC
    for c in ["provided_rank_score (max over alts)", "AVI phred (Atlas, max over alts)"]:
        sc = R[c].values[incat]
        jk = [fast_auc(np.delete(sc, i), np.delete(lab, i)) for i in range(len(sc))]
        add("familywise", c, "leave-one-variant-out AUC range (replicated label): median / min / max", float(np.median(jk)),
            float(np.min(jk)), float(np.max(jk)), np.nan, len(sc), int(n1), int(len(lab) - n1))
    # ---- top-k enrichment of prior-evidence positives among AG top 10 (provided rank score)
    from scipy.stats import hypergeom
    top = R["provided_rank_score (max over alts)"].sort_values(ascending=False).index[:10]
    for ln, lab in [("replicated >=2 pubs", t.lab_repl), ("GW-sig >=1 study", t.lab_gws), ("seed", t.seed)]:
        N, K, k = len(t), int(lab.sum()), int(lab.loc[top].sum())
        add("topk", "provided_rank_score top-10", f"count of '{ln}' in top 10 (expected {10*K/N:.1f})", k, np.nan, np.nan,
            hypergeom.sf(k - 1, N, K, 10), 10, K, N - K, note="hypergeometric upper-tail p")

    # ---- AG vs older predictors (agreement between methods, and paired AUC difference)
    agp = R["provided_rank_score (max over alts)"].values
    for nm, col in [("phyloP100way", t.phyloP100.values), ("phastCons100way", t.phastCons100.values)]:
        rho, lo, hi, p, n = spearman_ci(agp, col)
        add("agree", f"AG rank score vs {nm}", "Spearman across variants", rho, lo, hi, p, n)
    for nm, col in [("ENCODE cCRE regulatory-class overlap", t.ccre_regulatory), ("ENCODE cCRE any overlap", t.ccre_overlap)]:
        A, lo, hi, p, npos, nneg = auc_ci(agp, col.values.astype(bool))
        add("agree", f"AG rank score vs {nm}", "AUC of AG rank for cCRE-overlap variants", A, lo, hi, p, npos + nneg, npos, nneg)
    tp = lambda s: set(pd.Series(s, index=ids).dropna().sort_values(ascending=False).index[:10])
    t10 = {"AG": tp(agp), "phyloP": tp(t.phyloP100.values),
           "AVI": tp(R["AVI phred (Atlas, max over alts)"].values) if avi is not None else set()}
    for a1, a2 in [("AG", "phyloP"), ("AG", "AVI"), ("phyloP", "AVI")]:
        ov = len(t10[a1] & t10[a2])
        add("agree", f"top-10 overlap {a1} vs {a2}", "count of shared variants (expected ~2.0 by chance)", ov, np.nan, np.nan,
            hypergeom.sf(ov - 1, len(ids), 10, 10), 10)
    # paired bootstrap: AUC(AlphaGenome score) - AUC(older predictor), replicated / GW-sig labels, all and non-coding only
    avp = R["AVI phred (Atlas, max over alts)"].values if avi is not None else None
    olds = [("phyloP100way", t.phyloP100.values), ("phastCons100way", t.phastCons100.values),
            ("ENCODE cCRE overlap", t.ccre_overlap.values.astype(float))]
    for ln in ["replicated: GW-sig in >=2 publications (vs other in-catalog)", "GW-significant in >=1 study (vs in-catalog p>=5e-8)"]:
        lab, mask0 = labs[ln]
        lab = lab.values.astype(bool)
        for sub, mask in [("all", mask0), ("non-coding only", mask0 & ~t.coding_like.values)]:
            for an, ag_s in [("AG rank score", agp), ("AVI phred", avp)]:
                if ag_s is None:
                    continue
                if sub != "all":  # non-coding-only absolute AUCs
                    for nm, col in [(an, ag_s)] + olds:
                        A, lo, hi, p, npos, nneg = auc_ci(col[mask], lab[mask])
                        add("noncoding", f"{nm} [{sub}]", "AUC: " + ln, A, lo, hi, p, npos + nneg, npos, nneg)
                for nm, col in olds:
                    s1, s2, l = ag_s[mask], col[mask], lab[mask]
                    d0 = fast_auc(s1, l) - fast_auc(s2, l)
                    lo, hi = boot(lambda i: fast_auc(s1[i], l[i]) - fast_auc(s2[i], l[i]), len(l))
                    add("agree", f"AUC({an}) - AUC({nm}) [{sub}]", ln, d0, lo, hi, np.nan, int(mask.sum()))

    out = pd.DataFrame(rows).drop_duplicates(['section', 'predictor_or_recipe', 'outcome'])
    out.to_csv(ROOT / "results/alphagenome/benchmark_vs_prior.csv", index=False)
    var = t.join(R.add_prefix("AG::"))
    if avi is not None:
        var["AVI_phred_max"] = avi.groupby("rs_id").avi_phred.max().reindex(var.index)
    var.to_csv(ROOT / "results/alphagenome/benchmark_vs_prior_variants.csv")
    if len(cod):
        cod.to_csv(ROOT / "results/alphagenome/benchmark_coding_predictors.csv")
    write_report(out, t, R, avi, cod, best, worst, fw_p, older_ok)


VERDICT = """## Plain verdict (written against the seeded run; re-check if inputs change)

- **Ad hoc rank score: weak, statistically unconvincing concordance.** Spearman with -log10 p is 0.12 [-0.20, 0.44]. It does lean the right way for replicated loci (AUC 0.65 [0.39, 0.89]; 4 of its top 10 are replicated lead SNPs vs 1.8 expected, p=0.065) but every interval includes 0.5. Across 23 recipe variants (per-track, per-cell-type, rank-average vs max-of-ranks, chromatin vs expression) Spearman ranges only 0.01 to 0.21 and every CI includes 0; replicated-label AUC ranges 0.57 to 0.88 depending on recipe, so the headline depends on the recipe and 44 variants cannot pick one.
- **GW-significant AUC 0.80 is not evidence of success:** only 3 of 44 catalog variants are non-GW-significant, so that AUC has a CI of 0.53 to 1.00 and p=0.097. Single recipes reach 0.89 to 0.98 (ATAC-only, myofibroblast-only) but those are the best of many looks; do not quote them.
- **Atlas AVI score (obtained free via the documented Atlas API, 70/70 variant-alleles) is the one recipe with a signal:** AUC 0.88 [0.70, 1.00] for replicated lead SNPs (n 9 vs 35), family-wise permutation p about 0.001 over all 23 recipes, stable under leave-one-out (0.87 to 0.95) and when the 3 missense variants are removed (0.85 [0.62, 1.00]). It does not track p-value size (Spearman 0.21 [-0.11, 0.50]), GW-significance (0.59), seeds (0.45) or GTEx eQTL (0.44). This was tested after seeing the other recipes, so treat it as hypothesis-generating until repeated on new data.
- **Older, free predictors do about as well on the same label:** ENCODE cCRE overlap AUC 0.71 [0.53, 0.89], phastCons 0.77 [0.54, 1.00], phyloP 0.68 [0.38, 0.98]. AVI beats them by 0.11 to 0.20 AUC with CIs that barely exclude 0 (non-coding only: +0.13 [-0.01, 0.31] vs cCRE). AVI's top 10 shares 7 of 10 variants with phyloP's top 10 (AVI contains conservation and AlphaMissense features), so its edge is not clearly independent of conservation. The ad hoc rank score is statistically indistinguishable from the older predictors (differences -0.03 to -0.12, CIs span 0).
- **Agreement between methods:** the rank score tracks ENCODE cCRE overlap very strongly (AUC 0.97 [0.92, 1.00]) but this is expected because AlphaGenome is trained on ENCODE-type data, so it is consistency with a regulatory annotation, not independent validation. It is uncorrelated with phyloP (rho 0.03) and with AVI at variant level only partly (4 of top 10 shared). Coding variants disagree: PHLDA3 rs35383942 is AG rank 4 but SIFT deleterious / PolyPhen benign / AlphaMissense likely-benign (0.20); PHLDA3 rs192314256 is AG rank 22 with SIFT 0.00, PolyPhen 0.68, AlphaMissense 0.91.
- **Net:** the AlphaGenome rank score is only weakly and non-significantly concordant with prior GWAS evidence and is no better than simple free annotations. The Atlas AVI score shows a replicated, multiplicity-corrected but small-sample association with replication status of the lead SNPs, which is encouraging but not proof of causal-variant detection. Neither supports claims of accuracy; concordance is modest at best."""

LIMITS = """## Limitations

- **Small n and few negatives:** 44 catalog variants, 9 replicated, 3 non-GW-significant. CIs are wide, bootstrap CIs on binary or tied predictors (cCRE, 3 negatives) are unreliable, and 7 labels x 23 recipes were examined; only the AUC and Spearman maxima were corrected for multiplicity (permutation), nothing else.
- **No ground truth:** none of these labels marks a causal variant. GWAS p-value and replication reflect power, allele frequency and LD as much as biology. "Replicated" counts publications, which may share biobank participants (Sakaue/Ishigaki both BBJ-based; Verma/Greene/Dand overlapping European-ancestry biobanks), so independence is partial.
- **LD:** lead SNPs are tags, not necessarily causal. A strong GWAS signal with a weak predicted effect (e.g. rs10863683, p=2e-79, low rank) is not a model failure; conversely a high score at a weak signal is not a true positive.
- **Rank score unvalidated:** the ad hoc rank-average has no calibration; it was built on the same variants and is not independent of choices made while looking at the data. Scores are variant-level maxima over alleles because the effect allele is the reference base for many variants (alt ambiguous), and over 6 skin conditions, which favours variants with any large single-condition effect.
- **Labels:** GTEx matches use only the mapped gene(s) (38 of 50 assessable), GTEx significance threshold is lenient (14 of 38 positive), and the literature-count label is open-access only and overlaps the seed list. Seed variants without a catalog p-value (6) are excluded from p-based tests. Zhu 2013 verified as Europe PMC PMID 23667473 (Chinese Han keloid association study); Nakashima 2010 verified through GWAS Catalog GCST000773 (PMID 20711176), not open access in Europe PMC.
- **Tools and access:** AVI scores were obtained with the documented Atlas client (`atlas.create`, `query_variant`, the project's alphagenome_variant_impact_score route) using the existing key; no scraping. The ENCODE SCREEN GraphQL coordinate search returned empty for a known-positive region, so cCRE overlap came from the ENCODE registry track served by UCSC. The GTEx by-location endpoint returned nothing for windows known to contain eQTLs, so per-gene eQTL lists were used. CADD was not run (no approved free route in the skill wrappers; calling its web API is a new service). AlphaMissense is shown for the 3 missense variants but is a 2023 deep-learning predictor, not an "older" one. UCSC phyloP/phastCons used 100-vertebrate alignment only.
- All results are predictions and descriptive statistics for research use only; no clinical or causal claims."""


def fmt(r):
    ci = "" if np.isnan(r.ci_lo) else f" [{r.ci_lo:.2f}, {r.ci_hi:.2f}]"
    p = "" if np.isnan(r.p_value) else (f", p={r.p_value:.2g}")
    return f"{r.estimate:.2f}{ci}{p}"


def write_report(out, t, R, avi, cod, best, worst, fw_p, older_ok):
    L = []
    w = L.append
    w("# Concordance of AlphaGenome 2 keloid-variant scores with prior evidence\n")
    w("Research only. There is no gold-standard causal-variant set for keloid, so this measures **concordance with "
      "prior, independent evidence**, not accuracy or success rate. Generated by `scripts/benchmark_vs_prior.py`; "
      "per-statistic numbers are in `results/alphagenome/benchmark_vs_prior.csv`, per-variant values in "
      "`benchmark_vs_prior_variants.csv`.\n")
    w(VERDICT + "\n")
    w("## Method\n")
    w(f"- **Variants:** {len(t)} (variant-level; scores are the max over alleles scored, 6 skin conditions and tracks). "
      f"{int(t.in_catalog.sum())} are in the GWAS Catalog keloid set; {int((~t.in_catalog).sum())} are literature/seed variants with no catalog p-value.")
    w("- **Prior-evidence labels (AlphaGenome never used):** (1) genome-wide significant (min catalog p < 5e-8); "
      "(2) replicated = GW-significant in >= 2 distinct publications (PMIDs; biobank overlap between studies is possible); "
      "(3) replicated in >= 3 publications; (4) literature-seed flag, and the subset Nakashima 2010 (PMID 20711176) / Zhu 2013; "
      "(5) GTEx v10 significant eQTL at the exact SNP, for the variant's mapped gene(s), in skin (sun exposed, not exposed) or cultured fibroblasts (38 of 50 variants assessable; the skill's by-location endpoint returned nothing for windows known to contain NEDD4 eQTLs, so per-gene eQTL lists were matched by position); "
      "(6) >= 3 Europe PMC open-access papers co-mentioning the rsID and keloid (open access only; overlaps with the seed list, so not independent of it). "
      "Continuous: -log10 of the smallest catalog p.")
    w("- **Tests:** Spearman with percentile-bootstrap 95% CI (5000 resamples of variants); AUC (= Mann-Whitney U / n1 n2) with bootstrap CI; "
      "top-10 hypergeometric enrichment. **Recipe sensitivity:** the ad hoc rank score (max and mean over alleles), "
      "a recomputed rank-average, a max-of-ranks, 8 single-track scores, 6 per-cell-type rankings, chromatin-only and expression-only averages, "
      "and the AlphaGenome Atlas AVI score when available. A family-wise permutation test (shuffle p across variants, take the max Spearman over all recipes) "
      "quantifies how much of the best recipe is selection luck (a second one does the same for the replicated-label AUC). RNA scores missing for a variant (no gene mask) were set to 0 (no predicted effect). Variant-level scores take the max over all alleles scored, because the effect allele is the reference base for many variants and the alt is therefore ambiguous.")
    w("- **Older comparators run free:** phyloP100way / phastCons100way at the SNP (UCSC); ENCODE cCRE overlap (ENCODE registry track via the UCSC API; "
      "the SCREEN GraphQL coordinate search returned empty even for a known positive HBB region so was not used); Ensembl VEP SIFT/PolyPhen/AlphaMissense for the missense variants. "
      "CADD was not run (not in the Ensembl skill wrapper; its own web API would be a new service, which needs approval).\n")

    def sect(sec, cols=("predictor_or_recipe", "outcome")):
        d = out[out.section == sec]
        w("| " + " | ".join(["Predictor / recipe", "Outcome", "n (pos/neg)", "Estimate [95% CI], p"]) + " |")
        w("|---|---|---|---|")
        for r in d.itertuples():
            npn = f"{int(r.n)}" + (f" ({int(r.n_pos)}/{int(r.n_neg)})" if r.n_pos == r.n_pos else "")
            w(f"| {r.predictor_or_recipe} | {r.outcome} | {npn} | {fmt(r)} |")
        w("")

    def pivot(sec):
        short = {"Spearman": "Spearman vs -log10 p", "AUC: GW-significant": "AUC GW-sig", "AUC: replicated: ": "AUC replicated>=2",
                 "AUC: literature seed": "AUC seed", "AUC: classic": "AUC Nakashima/Zhu seed", "AUC: GTEx": "AUC GTEx eQTL",
                 "AUC: EuropePMC": "AUC lit>=3 papers"}
        d = out[out.section == sec].copy()
        d["col"] = d.outcome.map(lambda o: next((v for k, v in short.items() if o.startswith(k)), o))
        cols = list(dict.fromkeys(d.col))
        w("| Predictor / recipe | " + " | ".join(cols) + " |")
        w("|---|" + "---|" * len(cols))
        for nm in dict.fromkeys(d.predictor_or_recipe):
            cells = []
            for c in cols:
                r = d[(d.predictor_or_recipe == nm) & (d.col == c)]
                cells.append(fmt(r.iloc[0]) if len(r) else "-")
            w(f"| {nm} | " + " | ".join(cells) + " |")
        w("")

    w("## Prior-evidence label counts\n")
    w("| Label | Positives | of |")
    w("|---|---|---|")
    w(f"| GW-significant (>=1 study) | {int(t.lab_gws.sum())} | {int(t.in_catalog.sum())} in catalog |")
    w(f"| Replicated >= 2 publications | {int(t.lab_repl.sum())} | {int(t.in_catalog.sum())} in catalog |")
    w(f"| Replicated >= 3 publications | {int(t.lab_repl3.sum())} | {int(t.in_catalog.sum())} in catalog |")
    w(f"| Literature seed | {int(t.seed.sum())} | {len(t)} |")
    w(f"| Seed Nakashima 2010 / Zhu 2013 | {int(t.seed_classic.sum())} | {len(t)} |")
    w(f"| GTEx skin/fibroblast eQTL at SNP (mapped gene) | {int(t.lab_gtex.fillna(False).sum())} | {int(t.gtex_assessed.sum())} assessable |")
    w(f"| Europe PMC >=3 OA papers co-mentioning rsID and keloid | {int(t.lab_lit.sum())} | {len(t)} |\n")
    w("## AlphaGenome recipes vs prior evidence\n")
    pivot("alphagenome")
    w("### Family-wise check across recipes\n")
    sect("familywise")
    w("## Non-coding variants only (3 missense variants removed)\n")
    pivot("noncoding")
    w("## Top 10 variants by method\n")
    cand = {"AG rank score": R["provided_rank_score (max over alts)"],
            "AVI phred": R["AVI phred (Atlas, max over alts)"] if avi is not None else None, "phyloP100": t.phyloP100}
    w("| Method | Top 10: rsID (publications GW-sig / min p) |")
    w("|---|---|")
    for k, v in cand.items():
        if v is None:
            continue
        top10 = v.sort_values(ascending=False).index[:10]
        w(f"| {k} | " + "; ".join(f"{i} ({int(t.loc[i, 'n_pmids_gws'])} / {t.loc[i, 'min_p']:.0e})" for i in top10) + " |")
    w("")
    w("## Top-10 enrichment (ad hoc rank score)\n")
    sect("topk")
    w("## Older predictors vs prior evidence\n")
    pivot("older")
    w("## Agreement between AlphaGenome and older predictors\n")
    sect("agree")
    if len(cod):
        w("## Coding variants: older predictors vs AlphaGenome\n")
        rk = R["provided_rank_score (max over alts)"].rank(ascending=False)
        w("| Variant | SIFT (min score, call) | PolyPhen (max, call) | AlphaMissense (max) | AG rank (of 50) | GWAS min p |")
        w("|---|---|---|---|---|---|")
        for rs, r in cod.iterrows():
            w(f"| {rs} | {r.get('sift_min', np.nan):.2f} {r.get('sift_call', '')} | {r.get('polyphen_max', np.nan):.2f} {r.get('polyphen_call', '')} | "
              f"{r.get('alphamissense_max', np.nan):.2f} | {int(rk[rs])} | {t.loc[rs, 'min_p']:.1e} |")
        w("")
    w(LIMITS)
    open(ROOT / "reports/benchmark_vs_prior.md", "w", encoding="utf-8").write("\n".join(L))


if __name__ == "__main__":
    main()
