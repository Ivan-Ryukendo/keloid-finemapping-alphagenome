"""MYO7A WT vs V679I AF3 comparison. Same method as compare_structures.py (copied, not imported: that module runs PHLDA3 at import)."""
import json, itertools, csv
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

R = Path(__file__).resolve().parent.parent / "results" / "alphafold"
WT = ("myo7a_wt_q13402", "fold_myo7a_wt_q13402")
VAR = ("myo7a_v679i_rs35641839", "fold_myo7a_v679i_rs35641839")
POS, AA_W, AA_V, HW = 679, "V", "I", 25
MOTOR = (65, 741)  # UniProt Q13402 "Myosin motor" domain (verified via UniProt + InterPro IPR036106/PS51456)
P = MMCIFParser(QUIET=True)

def load(folder, prefix, i):
    d = R / folder
    s = P.get_structure("m", str(d / f"{prefix}_model_{i}.cif"))
    atoms = list(s.get_atoms())
    full = json.load(open(d / f"{prefix}_full_data_{i}.json"))
    assert len(atoms) == len(full["atom_plddts"]) and len(list(s.get_chains())) == 1
    ca, plddt, seq = [], [], ""
    for k, a in enumerate(atoms):
        if a.get_id() == "CA":
            ca.append(a.coord); plddt.append(full["atom_plddts"][k])
            seq += index_to_one(three_to_index(a.get_parent().get_resname()))
    rank = json.load(open(d / f"{prefix}_summary_confidences_{i}.json"))["ranking_score"]
    return dict(ca=np.array(ca), plddt=np.array(plddt), seq=seq, pae=np.array(full["pae"]), rank=rank)

def fit(mob, ref, idx=slice(None)):
    """Kabsch fit on idx only; return (rotated mob over all residues, RMSD over idx)."""
    mc, rc = mob[idx].mean(0), ref[idx].mean(0)
    u, _, vt = np.linalg.svd((mob[idx] - mc).T @ (ref[idx] - rc))
    rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    sup = (mob - mc) @ rot + rc
    return sup, float(np.sqrt(((sup[idx] - ref[idx]) ** 2).sum(1).mean()))

W = slice(POS - HW - 1, POS + HW)          # residues 654-704
MS = slice(MOTOR[0] - 1, MOTOR[1])          # motor domain
W5 = slice(POS - 6, POS + 5)                # 674-684, as PHLDA3 +-5

def pair(a, b):
    sup, glob = fit(a["ca"], b["ca"])
    wg = float(np.sqrt(((sup[W5] - b["ca"][W5]) ** 2).sum(1).mean()))  # +-5 under global fit (PHLDA3 metric)
    return dict(global_ca=glob, win5_globalfit=wg, win25_refit=fit(a["ca"], b["ca"], W)[1], motor_refit=fit(a["ca"], b["ca"], MS)[1])

METRICS = [("global_ca", "Whole-chain CA RMSD"), ("win5_globalfit", f"CA RMSD {POS-5}-{POS+5}, global fit"),
           ("win25_refit", f"CA RMSD {POS-HW}-{POS+HW}, refit on window"), ("motor_refit", f"CA RMSD motor domain {MOTOR[0]}-{MOTOR[1]}, refit on domain")]
rng = lambda xs: (float(np.mean(xs)), float(min(xs)), float(max(xs)))
fmt = lambda t: f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})"

wt = [load(*WT, i) for i in range(5)]
v = [load(*VAR, i) for i in range(5)]
best = max(range(5), key=lambda i: wt[i]["rank"])
assert len(wt[0]["seq"]) == 2215 and all(len(m["seq"]) == 2215 for m in wt + v)
assert all(m["seq"][POS - 1] == AA_W for m in wt) and all(m["seq"][POS - 1] == AA_V for m in v)
assert all(m["seq"][:POS-1] + AA_W + m["seq"][POS:] == wt[0]["seq"] for m in v), "other differences"
assert all(m["seq"] == wt[0]["seq"] for m in wt)

var = [pair(m, wt[best]) for m in v]
floor = [pair(wt[i], wt[j]) for i, j in itertools.combinations(range(5), 2)]
rows, md_rows = [], []
for key, label in METRICS:
    vv, ff = rng([p[key] for p in var]), rng([p[key] for p in floor])
    det = vv[0] > ff[2]
    rows.append(dict(metric=label, variant_mean=vv[0], variant_min=vv[1], variant_max=vv[2], wtwt_mean=ff[0], wtwt_min=ff[1], wtwt_max=ff[2], exceeds_wtwt_max=det))
    md_rows.append(f"| {label} (A) | {fmt(vv)} | {fmt(ff)} | {'yes' if det else 'no'} |")

def info(ms):
    return dict(mean_plddt=np.mean([m["plddt"].mean() for m in ms]), frac70=np.mean([(m["plddt"] < 70).mean() for m in ms]),
                mean_pae=np.mean([m["pae"].mean() for m in ms]),
                loc_pae=np.mean([(m["pae"][POS-1].mean() + m["pae"][:, POS-1].mean()) / 2 for m in ms]),
                p=np.mean([m["plddt"] for m in ms], 0), win_plddt=np.mean([m["plddt"][W].mean() for m in ms]))
iw, iv = info(wt), info(v)
per_model = lambda ms, tag: [f"| {tag} {i} | {m['rank']:.3f} | {m['plddt'].mean():.1f} | {(m['plddt']<70).mean()*100:.1f}% | {m['plddt'][POS-1]:.1f} |" for i, m in enumerate(ms)]
rows2 = [dict(metric=f"mean pLDDT (all {len(wt[0]['seq'])} res), WT|variant", variant_mean=iv["mean_plddt"], wtwt_mean=iw["mean_plddt"]),
         dict(metric="fraction pLDDT<70, WT|variant", variant_mean=iv["frac70"], wtwt_mean=iw["frac70"]),
         dict(metric=f"mean pLDDT window {POS-HW}-{POS+HW}, WT|variant", variant_mean=iv["win_plddt"], wtwt_mean=iw["win_plddt"]),
         dict(metric="mean PAE (A), WT|variant", variant_mean=iv["mean_pae"], wtwt_mean=iw["mean_pae"]),
         dict(metric=f"local PAE at {POS} (row+col mean, A), WT|variant", variant_mean=iv["loc_pae"], wtwt_mean=iw["loc_pae"])]
rr = range(POS - 4, POS + 3)
for r in rr:
    rows2.append(dict(metric=f"pLDDT residue {r+1}, WT|variant", variant_mean=iv["p"][r], wtwt_mean=iw["p"][r]))
cols = ["metric", "variant_mean", "variant_min", "variant_max", "wtwt_mean", "wtwt_min", "wtwt_max", "exceeds_wtwt_max"]
with open(R / "structure_comparison_myo7a.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, restval=""); w.writeheader()
    for r in rows + rows2:
        w.writerow({k: (round(x, 3) if isinstance(x, float) else x) for k, x in r.items()})
# in rows2, "wtwt_mean" column holds the WT value and "variant_mean" the variant value (labelled in metric)

any_det = any(r["exceeds_wtwt_max"] for r in rows)
verdict = ("at least one metric exceeds the WT-vs-WT maximum (see table); still not evidence of a functional effect" if any_det
           else "no detectable difference: all variant means are within the WT-vs-WT seed noise")
md = f"""# MYO7A AlphaFold 3 structure comparison: WT vs V679I (rs35641839)

AlphaFold Server predictions, 5 models per job, 2,215 aa (UniProt Q13402). Research only.
Reference: WT model {best} (top ranking_score). Residue check passed: V{POS} in all 5 WT models, I{POS} in all 5 variant models, otherwise identical sequence.

**Domain of residue {POS}:** Myosin motor domain (UniProt feature 65-741; InterPro IPR036106 79-729, PS51456 65-741, Pfam PF00063 67-729), in the C-terminal part of the motor, before the first IQ motif (745-765). It lies outside the UniProt actin-binding region (632-639). CATH places it in 6.20.240.20 (671-743).

| Metric | Variant vs WT-best, mean (range) | WT-vs-WT noise, 10 pairs, mean (range) | Variant mean > WT-vs-WT max |
|---|---|---|---|
""" + "\n".join(md_rows) + f"""

Whole-chain RMSD is dominated by domain/tail movement and is a poor local readout; use the refit rows.

| Confidence | WT (mean of 5) | Variant (mean of 5) |
|---|---|---|
| Mean pLDDT, all residues | {iw['mean_plddt']:.1f} | {iv['mean_plddt']:.1f} |
| Fraction of residues pLDDT < 70 | {iw['frac70']*100:.1f}% | {iv['frac70']*100:.1f}% |
| Mean pLDDT, residues {POS-HW}-{POS+HW} | {iw['win_plddt']:.1f} | {iv['win_plddt']:.1f} |
| Mean PAE, whole chain (A) | {iw['mean_pae']:.2f} | {iv['mean_pae']:.2f} |
| Local PAE at residue {POS} (row+column mean, A) | {iw['loc_pae']:.2f} | {iv['loc_pae']:.2f} |

| Residue | {' | '.join(str(r+1) for r in rr)} |
|---|{'---|'*7}
| pLDDT WT | {' | '.join(f'{iw["p"][r]:.1f}' for r in rr)} |
| pLDDT variant | {' | '.join(f'{iv["p"][r]:.1f}' for r in rr)} |

Per model (rank score, mean pLDDT, % residues <70, pLDDT at {POS}):

| Model | ranking_score | mean pLDDT | % <70 | pLDDT {POS} |
|---|---|---|---|---|
""" + "\n".join(per_model(wt, "WT") + per_model(v, "Variant")) + f"""

**Verdict: {verdict}.** Rule: variant mean must exceed the maximum WT-vs-WT value for the same metric.

## Caveats

- AlphaFold 3 does not predict the stability or binding effect of point mutations; small RMSD/pLDDT shifts are not evidence of functional impact, and a null result does not show the variant is harmless.
- V679I is a conservative change (Val to Ile) predicted benign/tolerated by VEP.
- n=5 models per condition, single run each; the noise floor is a rough estimate. This is a single-chain prediction without calmodulin/light chains, actin or nucleotide.
"""
(R / "structure_comparison_myo7a.md").write_text(md, encoding="utf-8")
print(md)
