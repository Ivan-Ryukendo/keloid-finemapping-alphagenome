"""Make README/paper figures (PNG 200 dpi + SVG). Run from repo root or anywhere:
python scripts/make_figures.py [--out DIR]"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=str(ROOT / "figures"))
OUT = Path(ap.parse_args().out)
OUT.mkdir(parents=True, exist_ok=True)
AG = ROOT / "results" / "alphagenome"

# Mid-tone colours with a white figure background: legible on light and dark README themes.
BLUE, GREY, INK, MUTED = "#2B6CB0", "#9AA5B1", "#1F2933", "#52606D"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white"})

def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)

top = pd.read_csv(AG / "finemap_top_per_locus.csv")
ranked = top[top.status == "ok"].copy()

# ---- Fig 1: workflow diagram
steps = [
    ("GWAS Catalog", "keloid trait MONDO_0005348\n39 lead loci"),
    ("LD fine-mapping", "39 leads -> 349 candidates\n(278 unique SNVs)"),
    ("AlphaGenome 2 scoring", "6 skin cell types\n1,644 / 1,668 pairs scored\n(4 variants failed: patch contigs)"),
    ("Rank-average", "7 regulatory metrics\n-> rank_score"),
    ("Per-locus top variant", "35 loci ranked"),
    ("GSE163973 scRNA", "mapped-gene expression\nby cell type"),
    ("Novelty check", "GWAS Catalog, GTEx,\nClinVar, literature"),
    ("Lab validation", "hypotheses for\nwet-lab collaborators"),
]
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis("off"); ax.set_xlim(0, 4); ax.set_ylim(0, 2.3)
pos = [(i % 4, 1.5 if i < 4 else 0.35) for i in range(8)]
if True:  # row 2 runs right-to-left so arrows follow the snake
    pos = [(i, 1.5) for i in range(4)] + [(3 - i, 0.35) for i in range(4)]
for i, ((t, s), (x, y)) in enumerate(zip(steps, pos)):
    last = i == 7
    ax.add_patch(FancyBboxPatch((x + .06, y), .88, .65, boxstyle="round,pad=0.01,rounding_size=0.04",
                                fc="#EBF4FF" if not last else "#FFF4E0", ec=BLUE if not last else "#B7791F", lw=1.4))
    ax.text(x + .5, y + .5, f"{i+1}. {t}", ha="center", va="center", fontweight="bold", fontsize=10)
    ax.text(x + .5, y + .22, s, ha="center", va="center", fontsize=8, color=MUTED)
    if i < 7:
        nx, ny = pos[i + 1]
        if ny == y:
            d = 1 if nx > x else -1
            ax.annotate("", xy=(nx + .5 - d * .44, ny + .325), xytext=(x + .5 + d * .44, y + .325),
                        arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
        else:
            ax.annotate("", xy=(nx + .5, ny + .65), xytext=(x + .5, y), arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
ax.set_title("Keloid variant-prioritisation pipeline (computational predictions; hypothesis-generating)", loc="left", fontsize=11)
save(fig, "fig1_pipeline_workflow")

# ---- Fig 2: top-10 loci
t10 = ranked.sort_values("rank_score", ascending=False).head(10).iloc[::-1]
labels = [f"{r.lead_rsid}  ->  {r.top_rsid}" for r in t10.itertuples()]
cols = [BLUE if s else GREY for s in t10.skin_accessible.astype(bool)]
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.barh(labels, t10.rank_score, color=cols, height=0.6)
for y, v in enumerate(t10.rank_score):
    ax.text(v + .005, y, f"{v:.3f}", va="center", fontsize=8, color=MUTED)
ax.set_xlim(0, 1.1); ax.set_xlabel("rank_score (rank-average of 7 AlphaGenome metrics; higher = stronger predicted effect)")
ax.set_title("Top 10 loci: best candidate variant per locus (lead -> top SNV)", loc="left", fontsize=11)
ax.legend(handles=[Line2D([0], [0], marker="s", ls="", color=BLUE, label="in skin-accessible chromatin"),
                   Line2D([0], [0], marker="s", ls="", color=GREY, label="not skin-accessible")],
          loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
save(fig, "fig2_top10_loci")

# ---- Fig 3: expression heatmap (GWAS-mapped genes, all conditions pooled)
g = pd.read_csv(AG / "gse163973_gene_by_celltype.csv")
g = g[(g.role == "gwas_mapped") & (g.condition == "all") & (g.cell_type != "Unknown")]
m = g.pivot(index="gene", columns="cell_type", values="mean_log1p_cp10k")
m = m.loc[m.max(axis=1).sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(m.values, aspect="auto", cmap="Blues")
ax.set_xticks(range(m.shape[1]), m.columns, rotation=40, ha="right")
ax.set_yticks(range(m.shape[0]), m.index, fontsize=8)
ax.spines[:].set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=.6); cb.set_label("mean log1p(CP10K)")
ax.set_title("GWAS-mapped genes: expression by cell type (GSE163973, all samples)", loc="left", fontsize=11)
save(fig, "fig3_gene_expression_heatmap")

# ---- Fig 4: score distribution + coverage
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1.6, 1]})
a.hist(ranked.rank_score, bins=12, color=BLUE, edgecolor="white")
a.set_xlabel("rank_score of per-locus top variant"); a.set_ylabel("loci")
a.set_title(f"Score distribution ({len(ranked)} ranked loci)", loc="left", fontsize=11)
cov = top.status.value_counts()
names = {"ok": "ranked", "no candidates (multiallelic)": "no candidates\n(multiallelic)", "candidates not scored": "candidates\nnot scored"}
b.bar([names[k] for k in cov.index], cov.values, color=[BLUE, GREY, GREY], width=.55)
for i, v in enumerate(cov.values):
    b.text(i, v + .5, str(v), ha="center")
b.set_ylabel("loci"); b.set_title(f"Locus coverage ({len(top)} leads)", loc="left", fontsize=11)
fig.tight_layout()
save(fig, "fig4_score_distribution_coverage")
# ---- Fig 5: GTEx skin eQTL (exact variant) for top variant per lead, colour = variant-level tier
import numpy as np
nv = pd.read_csv(ROOT / "results" / "novelty" / "novelty_all_loci.csv")
e = nv[nv.gtex_skin_eqtl == "y"].copy()
e["nlp"] = -np.log10(e.gtex_skin_p); e = e.sort_values("nlp")
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.barh([f"{r.lead} -> {r.variant}" for r in e.itertuples()], e.nlp,
        color=[BLUE if t == 2 else GREY for t in e.tier], height=0.6)
for y, (v, p) in enumerate(zip(e.nlp, e.gtex_skin_p)):
    ax.text(v + .3, y, f"p={p:.1e}", va="center", fontsize=8, color=MUTED)
ax.set_xlim(0, e.nlp.max() * 1.2); ax.set_xlabel("-log10 p, best GTEx skin / fibroblast eQTL for the exact variant")
ax.set_title("Top variants with a GTEx skin eQTL (blue = tier 2, grey = tier 1)", loc="left", fontsize=11)
save(fig, "fig5_gtex_skin_eqtl")
print("wrote", sorted(p.name for p in OUT.iterdir()))
