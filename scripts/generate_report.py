"""Assemble reports/keloid_findings_summary.md from existing result files (no new analysis).
Sections that depend on a missing file are marked PENDING, never silently dropped."""
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "results"
rank = pd.read_csv(R / "alphagenome" / "variant_ranking.csv")
gwas = pd.read_csv(ROOT / "data" / "raw" / "keloid_gwas_variants.csv").drop_duplicates("rs_id")
top = rank.head(10).merge(gwas[["rs_id", "p_value", "in_gwas_catalog_keloid"]], on="rs_id", how="left")


def md_table(df):
    return df.to_markdown(index=False, floatfmt=".3g")


def part(path, title):
    p = ROOT / path
    return p.read_text(encoding="utf-8") if p.exists() else f"_PENDING: {title} ({path} not found)_"


struct = []
for f in ["structure_comparison.csv", "structure_comparison_myo7a.csv"]:
    p = R / "alphafold" / f
    if p.exists():
        d = pd.read_csv(p)
        struct.append(d[[c for c in ["variant", "ca_rmsd_mean", "noise_ca_mean", "local_rmsd_mean", "noise_local_mean",
                                     "verdict"] if c in d.columns]])
struct_md = md_table(pd.concat(struct)) if struct else "_PENDING: structure comparison_"

out = f"""# Keloid computational findings (research only, not clinical)
Generated {date.today()}. AlphaGenome 2 (API) + AlphaFold 3 (AlphaFold Server). Predictions, not validation.

## 1. Regulatory scoring (AlphaGenome 2, 6 skin cell conditions)
50 GWAS/literature variants, 70 variant-allele pairs. Rank score = average of per-track percentile ranks; it is an
ad hoc summary, **not** a validated pathogenicity score. Top 10:

{md_table(top[["rs_id", "alt", "gene", "consequence", "rank_score", "p_value", "in_gwas_catalog_keloid"]])}

## 2. Structure (AlphaFold 3, wild-type vs variant)
A change counts only if it exceeds the wild-type seed-to-seed noise.

{struct_md}

## 3. Agreement with prior evidence
{part("reports/benchmark_vs_prior.md", "benchmark")}

## 4. Biological context
{part("reports/context_notes.md", "context notes")}

## 5. Limitations
- Most keloid GWAS hits are non-coding and lead SNPs tag nearby causal variants (LD).
- No gold-standard causal variant set exists; scores show concordance, not accuracy.
- AlphaFold does not predict stability or binding effects of single substitutions.
- AlphaFold Server outputs are non-commercial; see results/alphafold/terms_of_use.md.
"""
(ROOT / "reports" / "keloid_findings_summary.md").write_text(out, encoding="utf-8")
assert "## 5." in out
print("wrote reports/keloid_findings_summary.md,", out.count("PENDING"), "pending section(s)")
