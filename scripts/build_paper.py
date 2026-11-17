"""Build paper/keloid_finemapping_paper.md (fills tables from the CSVs) and the PDF.
python scripts/build_paper.py [--pdf]   (needs: markdown; Chrome or Edge for --pdf)
Tables are generated from results/ so quoted numbers cannot drift from the data."""
import argparse, base64, html, re, subprocess, sys, tempfile
from pathlib import Path
import markdown
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
P = ROOT / "paper"
ap = argparse.ArgumentParser(); ap.add_argument("--pdf", action="store_true"); ap.add_argument("--browser", default="")
args = ap.parse_args()

nv = pd.read_csv(ROOT / "results/novelty/novelty_all_loci.csv")
top = pd.read_csv(ROOT / "results/alphagenome/finemap_top_per_locus.csv").set_index("lead_rsid")
rk = pd.read_csv(ROOT / "results/alphagenome/finemap_ranking_all.csv")
gene = lambda g: str(g).split(" | ")[0].replace(";", "; ")
fp = lambda p: "" if pd.isna(p) else f"{p:.2e}".replace("e-0", "e-")

def eqtl(r):
    return {"y": f"y ({fp(r.gtex_skin_p)})", "n": "n", "not_in_gtex": "nt"}.get(r.gtex_skin_eqtl, str(r.gtex_skin_eqtl))

# Table 1: top ten per-locus variants
t1 = ["| # | Lead → variant | Gene (Catalog) | rank_score | max \\|effect\\| | Skin chromatin flag | r² EUR / EAS | GTEx skin eQTL | Tier (literal) |", "|---|---|---|---|---|---|---|---|---|"]
for i, r in enumerate(nv.dropna(subset=["rank_score"]).sort_values("rank_score", ascending=False).head(10).itertuples(), 1):
    t = top.loc[r.lead]
    x = rk[(rk.lead_rsid == r.lead) & (rk.rs_id == r.variant)].iloc[0]
    r2 = lambda v: "–" if pd.isna(v) else f"{v:.2f}"
    flag = ("yes: " + str(t.skin_acc_tracks)) if t.skin_accessible else "**no**"
    t1.append(f"| {i} | {r.lead} → {r.variant} | {gene(r.mapped_gene)} | {t.rank_score:.3f} | {t.max_abs_effect:.2f} | {flag} | {r2(x.r2_EUR)} / {r2(x.r2_EAS)} | {eqtl(r)} | {r.tier} ({r.tier_ld_rule_literal}) |")

# Table S1: all 39 leads
s1 = ["| Lead → variant | Gene | rank_score | Tier | Literal | Depth | GTEx skin eQTL |", "|---|---|---|---|---|---|---|"]
for r in nv.sort_values("rank_score", ascending=False, na_position="last").itertuples():
    rs = "" if pd.isna(r.rank_score) else f"{r.rank_score:.3f}"
    s1.append(f"| {r.lead} → {r.variant} | {gene(r.mapped_gene)} | {rs} | {r.tier} | {r.tier_ld_rule_literal} | {'D' if r.evidence_depth == 'deep' else 'L'} | {eqtl(r)} |")

src = (P / "paper_source.md").read_text(encoding="utf8")
md = src.replace("{{T1}}", "\n".join(t1)).replace("{{S1}}", "\n".join(s1))
assert "{{" not in md
(P / "keloid_finemapping_paper.md").write_text(md, encoding="utf8")
print("wrote paper md;", len(nv), "leads;", nv.tier.value_counts().to_dict())

def build_html():
    body = markdown.markdown(md, extensions=["tables", "sane_lists"])
    def emb(m):
        f = (P / m.group(1)).resolve()
        return f'src="data:image/png;base64,{base64.b64encode(f.read_bytes()).decode()}"'
    body = re.sub(r'src="(\.\./figures/[^"]+)"', emb, body)
    css = """@page{size:A4;margin:16mm 15mm}body{font:10pt/1.4 'Segoe UI',Arial,sans-serif;color:#1f2933;max-width:100%}
h1{font-size:17pt;line-height:1.2}h2{font-size:13pt;border-bottom:1px solid #ccd;padding-bottom:2px;margin-top:18px;break-after:avoid}h3{font-size:11pt;break-after:avoid}
img{max-width:100%;max-height:125mm;width:auto;display:block;margin:6px auto;break-inside:avoid}p>img+em,p{orphans:3;widows:3}
table{border-collapse:collapse;width:100%;font-size:7.8pt;margin:8px 0;break-inside:auto}th,td{border:1px solid #bbc;padding:2px 4px;vertical-align:top}th{background:#eef2f7}tr{break-inside:avoid}
blockquote{margin:8px 0;padding:4px 10px;border-left:3px solid #b7791f;background:#fff7e6}code{font-size:8.5pt}"""
    return f"<!doctype html><html><head><meta charset='utf-8'><title>Keloid fine-mapping preprint</title><style>{css}</style></head><body>{body}</body></html>"

(P / "keloid_finemapping_paper.html").write_text(build_html(), encoding="utf8")
if args.pdf:
    b = args.browser or next((c for c in ["C:/Program Files/Google/Chrome/Application/chrome.exe",
        "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe", "C:/Program Files/Microsoft/Edge/Application/msedge.exe"] if Path(c).exists()), "")
    assert b, "no Chrome/Edge found; pass --browser"
    prof = Path(tempfile.mkdtemp(prefix="chrome_prof_", dir=str(ROOT.parent / "Keloid_Research_Workspace" / "cache") if (ROOT.parent / "Keloid_Research_Workspace" / "cache").exists() else None))
    out = P / "keloid_finemapping_paper.pdf"
    subprocess.run([b, "--headless=new", "--disable-gpu", f"--user-data-dir={prof}", "--no-pdf-header-footer",
                    f"--print-to-pdf={out}", (P / "keloid_finemapping_paper.html").as_uri()], check=True, timeout=180)
    print("wrote", out, out.stat().st_size, "bytes")
