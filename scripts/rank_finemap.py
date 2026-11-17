"""Rank LD-proxy variants per keloid locus from finemap_scores_long.csv (read-only on inputs). Pre-registered, untuned.
Score = benchmark_vs_prior 'rank-average, 7 provided metrics': per metric take max over the 6 conditions
(RNA as |lfc|), percentile-rank across all scored variants, average. Atlas AVI not available for these variants
(not queried: new calls), so it is not used. Accessible flag: ATAC/DNASE/CHIP_HISTONE log2fc >= own 90th
percentile over ALL rows, in skin fibroblast or keratinocyte. Every lead is reported, incl. ones with no scored variant."""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
ROOT = Path(__file__).resolve().parent.parent
A, P = ROOT / "results/alphagenome", ROOT / "data/processed"
M7 = ["ATAC_log2fc_absmax", "CAGE_log2fc_absmax", "CHIP_HISTONE_log2fc_absmax", "CHIP_TF_log2fc_absmax",
      "DNASE_log2fc_absmax", "SC_ATAC_log2fc_absmax", "RNA_lfc_absmax"]
ACC = ["ATAC_log2fc_absmax", "DNASE_log2fc_absmax", "CHIP_HISTONE_log2fc_absmax"]
SKIN = ["UBERON:0002097_CL:0002620", "UBERON:0002097_CL:0000312"]  # skin fibroblast, skin keratinocyte

s = pd.read_csv(A / "finemap_scores_long.csv")
s["RNA_lfc_absmax"] = s.RNA_lfc_absmax.abs().fillna(0)
cand = pd.read_csv(P / "finemap_candidates.csv")
exc = pd.read_csv(P / "finemap_excluded.csv")
thr = {c: s[c].quantile(0.9) for c in ACC}
s["max_abs_effect"] = s[M7].max(axis=1)

key = ["rs_id", "alt"]
v = s.groupby(key)[M7].max()
v["rank_score"] = pd.DataFrame({c: stats.rankdata(v[c]) / len(v) for c in M7}, index=v.index).mean(axis=1)
v["max_abs_effect"] = v[M7].max(axis=1)
sk = s[s.condition_id.isin(SKIN)].groupby(key)
v["skin_accessible"] = pd.concat([sk[c].max() >= thr[c] for c in ACC], axis=1).any(axis=1).reindex(v.index).fillna(False)
v["skin_acc_tracks"] = pd.concat([(sk[c].max() >= thr[c]).map({True: c.split("_")[0] + ("_H" if "HISTONE" in c else ""), False: ""})
                                  for c in ACC], axis=1).apply(lambda r: "+".join(x for x in r if x), axis=1).reindex(v.index)
for cn, d in s.groupby("condition_name"):  # per-condition summary: same rank-average within condition
    x = d.set_index(key)[M7]
    v[f"rank_{cn}"] = pd.DataFrame({c: stats.rankdata(x[c]) / len(x) for c in M7}, index=x.index).mean(axis=1).reindex(v.index)
v = v.reset_index()

# expand to every lead that lists the variant (scoring run deduped by proxy+alt)
cc = cand.rename(columns={"proxy_rsid": "rs_id"})[["lead_rsid", "rs_id", "alt", "r2_EUR", "r2_EAS", "consequence"]]
al = cc.merge(v, on=key, how="inner")
al["is_lead"] = al.lead_rsid == al.rs_id
al["locus_rank"] = al.groupby("lead_rsid").rank_score.rank(ascending=False, method="min")
al = al.sort_values(["lead_rsid", "locus_rank"])
al.to_csv(A / "finemap_ranking_all.csv", index=False)

leads = sorted(set(cand.lead_rsid) | set(exc.lead_rsid))
top = al[al.locus_rank == 1].drop_duplicates("lead_rsid").set_index("lead_rsid")
n_cand = cand.groupby("lead_rsid").size(); n_sc = al.groupby("lead_rsid").size()
rows = []
for l in leads:
    r = dict(lead_rsid=l, n_candidates=int(n_cand.get(l, 0)), n_scored=int(n_sc.get(l, 0)),
             n_excluded=int((exc.lead_rsid == l).sum()))
    if l in top.index:
        t = top.loc[l]
        r.update(top_rsid=t.rs_id, top_alt=t.alt, top_is_lead=t.is_lead, r2_EUR=t.r2_EUR, consequence=t.consequence,
                 rank_score=t.rank_score, max_abs_effect=t.max_abs_effect, skin_accessible=t.skin_accessible,
                 skin_acc_tracks=t.skin_acc_tracks, status="ok")
    else:
        r["status"] = ("candidates not scored" if n_cand.get(l, 0) else "no candidates (" + ";".join(sorted(set(exc[exc.lead_rsid == l].reason))) + ")")
    rows.append(r)
tp = pd.DataFrame(rows).sort_values("rank_score", ascending=False, na_position="last")
tp.to_csv(A / "finemap_top_per_locus.csv", index=False)
print("thresholds (p90, all rows):", {k: round(x, 4) for k, x in thr.items()})
print("rows", len(s), "variants", len(v), "locus-variant pairs", len(al), "leads", len(leads), "with top", len(top))
assert len(top) == al.lead_rsid.nunique() and al.groupby("lead_rsid").locus_rank.min().eq(1).all()
print(tp.head(10)[["lead_rsid", "top_rsid", "rank_score", "skin_accessible", "n_scored"]].to_string())
print(tp[tp.status != "ok"][["lead_rsid", "status"]].to_string())
