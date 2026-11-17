"""Compare AF3 variant vs WT models (Calpha RMSD, local RMSD, pLDDT, PAE) with a WT-vs-WT seed noise floor."""
import json, itertools, csv
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

R = Path(__file__).resolve().parent.parent / "results" / "alphafold"
WT = ("phlda3_wt_q9y5j5", "fold_phlda3_wt_q9y5j5")
VARS = {  # name: (folder, prefix, position, expected WT aa, expected variant aa)
    "R28Q": ("phlda3_r28q_rs35383942", "fold_phlda3_r28q_rs35383942", 28, "R", "Q"),
    "E62G": ("phlda3_e62g_rs192314256", "fold_phlda3_e62g_rs192314256", 62, "E", "G"),
}
P = MMCIFParser(QUIET=True)

def load(folder, prefix, i):
    d = R / folder
    s = P.get_structure("m", str(d / f"{prefix}_model_{i}.cif"))
    atoms = list(s.get_atoms())
    plddt_all = json.load(open(d / f"{prefix}_full_data_{i}.json"))["atom_plddts"]
    assert len(atoms) == len(plddt_all)
    ca, plddt, seq = [], [], ""
    for k, a in enumerate(atoms):
        if a.get_id() == "CA":
            ca.append(a.coord); plddt.append(plddt_all[k])
            seq += index_to_one(three_to_index(a.get_parent().get_resname()))
    pae = np.array(json.load(open(d / f"{prefix}_full_data_{i}.json"))["pae"])
    rank = json.load(open(d / f"{prefix}_summary_confidences_{i}.json"))["ranking_score"]
    return dict(ca=np.array(ca), plddt=np.array(plddt), seq=seq, pae=pae, rank=rank)

def kabsch(mob, ref):
    """Return mob superposed onto ref (all points)."""
    mc, rc = mob.mean(0), ref.mean(0)
    u, _, vt = np.linalg.svd((mob - mc).T @ (ref - rc))
    rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    return (mob - mc) @ rot + rc

def rmsd(a, b, idx=slice(None)):
    return float(np.sqrt(((a[idx] - b[idx]) ** 2).sum(1).mean()))

def pair(a, b, pos):
    sup = kabsch(a["ca"], b["ca"])
    w = slice(max(pos - 6, 0), pos + 5)  # residues pos-5..pos+5 (1-based)
    return rmsd(sup, b["ca"]), rmsd(sup, b["ca"], w)

def stats(x):
    x = np.array(x); return f"{x.mean():.2f} ({x.min():.2f}-{x.max():.2f})"

wt = [load(*WT, i) for i in range(5)]
best = max(range(5), key=lambda i: wt[i]["rank"])
rows, md = [], []
for name, (fo, pre, pos, aa_w, aa_v) in VARS.items():
    v = [load(fo, pre, i) for i in range(5)]
    for m in wt: assert m["seq"][pos - 1] == aa_w, ("WT residue", name)
    for m in v: assert m["seq"][pos - 1] == aa_v, ("variant residue", name)
    assert all(m["seq"][:pos-1] + aa_w + m["seq"][pos:] == wt[0]["seq"] for m in v), "other differences"
    var = [pair(m, wt[best], pos) for m in v]
    floor = [pair(wt[i], wt[j], pos) for i, j in itertools.combinations(range(5), 2)]
    floor_best = [pair(wt[i], wt[best], pos) for i in range(5) if i != best]
    wp, vp = np.mean([m["plddt"] for m in wt], 0), np.mean([m["plddt"] for m in v], 0)
    rng = range(pos - 4, pos + 3)  # 0-based idx for pos-3..pos+3
    wpae, vpae = np.mean([m["pae"].mean() for m in wt]), np.mean([m["pae"].mean() for m in v])
    loc = lambda ms: np.mean([(m["pae"][pos-1].mean() + m["pae"][:, pos-1].mean()) / 2 for m in ms])
    plddt_rows = {f"plddt_res{r+1}": (round(wp[r], 1), round(vp[r], 1)) for r in rng}
    above = lambda x, f: float(np.mean(x)) > max(f)  # mean exceeds the whole noise range
    verdict = "detectable change" if above([a for a, _ in var], [a for a, _ in floor]) or above([b for _, b in var], [b for _, b in floor]) else "within seed noise"
    rows.append(dict(variant=name, wt_ref_model=best, ca_rmsd_mean=np.mean([a for a, _ in var]), ca_rmsd_min=min(a for a, _ in var), ca_rmsd_max=max(a for a, _ in var),
        local_rmsd_mean=np.mean([b for _, b in var]), local_rmsd_min=min(b for _, b in var), local_rmsd_max=max(b for _, b in var),
        noise_ca_mean=np.mean([a for a, _ in floor]), noise_ca_min=min(a for a, _ in floor), noise_ca_max=max(a for a, _ in floor),
        noise_local_mean=np.mean([b for _, b in floor]), noise_local_min=min(b for _, b in floor), noise_local_max=max(b for _, b in floor),
        mean_pae_wt=wpae, mean_pae_var=vpae, delta_pae=vpae - wpae, local_pae_wt=loc(wt), local_pae_var=loc(v), verdict=verdict,
        **{k: f"{a}|{b}" for k, (a, b) in plddt_rows.items()}))
    md.append(f"### {name} (residue {pos}: {aa_w}->{aa_v}, verified in all 5 WT and 5 variant CIFs)\n\n"
        f"Reference: WT model {best} (top ranking_score).\n\n| Metric | Variant vs WT-best, mean (range) | WT-vs-WT noise, 10 pairs, mean (range) |\n|---|---|---|\n"
        f"| Whole-chain CA RMSD (A) | {stats([a for a,_ in var])} | {stats([a for a,_ in floor])} |\n"
        f"| CA RMSD residues {pos-5}-{pos+5} (A) | {stats([b for _,b in var])} | {stats([b for _,b in floor])} |\n\n"
        f"WT-vs-WT noise vs the best WT model only (4 pairs): CA {stats([a for a,_ in floor_best])}, local {stats([b for _,b in floor_best])}.\n\n"
        f"Mean PAE: WT {wpae:.2f} -> variant {vpae:.2f} (delta {vpae-wpae:+.2f}). Local PAE at residue {pos} (row+column mean): WT {loc(wt):.2f} -> variant {loc(v):.2f}.\n\n"
        f"| Residue | " + " | ".join(str(r+1) for r in rng) + " |\n|---|" + "---|" * 7 + "\n"
        f"| pLDDT WT (mean of 5) | " + " | ".join(f"{wp[r]:.1f}" for r in rng) + " |\n"
        f"| pLDDT variant (mean of 5) | " + " | ".join(f"{vp[r]:.1f}" for r in rng) + " |\n\n"
        f"**Verdict: {verdict}.** (Rule: variant mean RMSD must exceed the maximum WT-vs-WT value.)\n")

with open(R / "structure_comparison.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]) + [k for r in rows for k in r if k not in rows[0]]); w.writeheader()
    for r in rows: w.writerow({k: (round(x, 3) if isinstance(x, float) else x) for k, x in r.items()})
(R / "structure_comparison.md").write_text(
    "# PHLDA3 AlphaFold 3 structure comparison (WT vs variants)\n\nAlphaFold Server predictions, 5 seeds/models per job. Research only. Seed-to-seed WT-vs-WT variation is the noise floor.\n\n"
    + "\n".join(md) + "\n## Caveats\n\n- AlphaFold 3 does not predict the stability or binding effect of point mutations; small RMSD/pLDDT shifts are not evidence of functional impact.\n"
    "- n=5 models per condition, single seed each; the noise floor is a rough estimate.\n- Both variants are missense changes in a mostly disordered/loop-rich small protein; interpret with wet-lab data.\n", encoding="utf-8")
print("\n".join(md))
