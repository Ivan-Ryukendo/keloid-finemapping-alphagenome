"""Fetch keloid GWAS associations (GWAS Catalog REST v2, trait MONDO_0005348) and flag
seed variants and likely coding variants. Output: data/raw/keloid_gwas_variants.csv"""
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
TRAIT = "MONDO_0005348"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "raw" / "keloid_gwas_variants.csv"
# Seed rsIDs from literature (see keloid_research_pipeline_claude_code_prompt.md section 4B)
SEEDS = {
    "rs873549": "Nakashima 2010 / Zhu 2013 (1q41)",
    "rs1442440": "Zhu 2013 (1q41)",
    "rs1511412": "Nakashima 2010 (3q22.3-23, FOXL2)",
    "rs940187": "Nakashima 2010 (3q22.3-23)",
    "rs8032158": "Nakashima 2010 (15q, NEDD4)",
    "rs2271289": "Zhu 2013 (NEDD4)",
    "rs747722": "African-American admixture study (MYO1E)",
    "rs28394564": "African-American admixture study (MYO1E)",
    "rs35641839": "African-American admixture study (MYO7A)",
}
CODING = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
          "inframe_insertion", "inframe_deletion", "synonymous_variant", "splice_donor_variant",
          "splice_acceptor_variant", "protein_altering_variant"}


def get(url, **params):
    for attempt in range(4):
        r = requests.get(url, params=params, timeout=60)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 404:
            return None
        time.sleep(2 * (attempt + 1))
    r.raise_for_status()


def fetch_associations():
    rows, page = [], 0
    while True:
        d = get(f"{BASE}/associations", efo_id=TRAIT, size=50, page=page)
        for a in d["_embedded"]["associations"]:
            alleles = a.get("snp_allele") or [{}]
            for al in alleles:
                rows.append({
                    "rs_id": al.get("rs_id"), "effect_allele": al.get("effect_allele"),
                    "p_value": a.get("p_value"), "beta_or_description": a.get("beta"),
                    "risk_frequency": a.get("risk_frequency"),
                    "location_grch38": ";".join(a.get("locations") or []),
                    "mapped_genes": ";".join(a.get("mapped_genes") or []),
                    "reported_trait": ";".join(a.get("reported_trait") or []),
                    "pubmed_id": a.get("pubmed_id"), "first_author": a.get("first_author"),
                    "study_accession": a.get("accession_id"), "association_id": a.get("association_id"),
                })
        if page + 1 >= d["page"]["totalPages"]:
            return rows
        page += 1


def annotate(rs_id):
    d = get(f"{BASE}/single-nucleotide-polymorphisms/{rs_id}")
    if not d:  # not in GWAS Catalog keloid set: fall back to Ensembl (GRCh38)
        e = get(f"https://rest.ensembl.org/variation/human/{rs_id}", **{"content-type": "application/json"})
        m = next((x for x in (e or {}).get("mappings", []) if x.get("assembly_name") == "GRCh38"), None)
        if not m:
            return {}
        return {"chrom": m["seq_region_name"], "pos_grch38": m["start"], "cytoband": None,
                "functional_class": None, "most_severe_consequence": e.get("most_severe_consequence"),
                "alleles": m.get("allele_string"), "maf": e.get("MAF")}
    loc = (d.get("locations") or [{}])[0]
    return {"chrom": loc.get("chromosome_name"), "pos_grch38": loc.get("chromosome_position"),
            "cytoband": (loc.get("region") or {}).get("name"),
            "functional_class": d.get("functional_class"),
            "most_severe_consequence": d.get("most_severe_consequence"),
            "alleles": d.get("alleles"), "maf": d.get("maf")}


def main():
    df = pd.DataFrame(fetch_associations())
    df = df.dropna(subset=["rs_id"])
    ids = sorted(set(df.rs_id) | set(SEEDS))
    rows = []
    for i in ids:
        rows.append({"rs_id": i, **annotate(i)})
        time.sleep(0.2)
    a = pd.DataFrame(rows)
    a["pos_grch38"] = a["pos_grch38"].astype("Int64")
    df = df.merge(a, on="rs_id", how="outer")
    df["in_gwas_catalog_keloid"] = df.association_id.notna()
    df["seed_source"] = df.rs_id.map(SEEDS)
    df["is_coding_like"] = df.most_severe_consequence.isin(CODING)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"associations: {int(df.in_gwas_catalog_keloid.sum())}, unique rsIDs: {df.rs_id.nunique()}")
    print("coding-like variants:", sorted(df[df.is_coding_like].rs_id.unique()))
    print("seed check:")
    for s, src in SEEDS.items():
        r = df[df.rs_id == s].iloc[0]
        print(f"  {s}: in_catalog={bool(r.in_gwas_catalog_keloid)} consequence={r.most_severe_consequence} "
              f"chr{r.chrom}:{r.pos_grch38}")
    assert len(df) > 0 and "pos_grch38" in df


if __name__ == "__main__":
    main()
