# PHLDA3 AlphaFold 3 structure comparison (WT vs variants)

AlphaFold Server predictions, 5 seeds/models per job. Research only. Seed-to-seed WT-vs-WT variation is the noise floor.

### R28Q (residue 28: R->Q, verified in all 5 WT and 5 variant CIFs)

Reference: WT model 0 (top ranking_score).

| Metric | Variant vs WT-best, mean (range) | WT-vs-WT noise, 10 pairs, mean (range) |
|---|---|---|
| Whole-chain CA RMSD (A) | 0.68 (0.44-0.83) | 0.82 (0.51-1.21) |
| CA RMSD residues 23-33 (A) | 0.34 (0.17-0.65) | 0.42 (0.19-0.75) |

WT-vs-WT noise vs the best WT model only (4 pairs): CA 0.69 (0.51-0.83), local 0.38 (0.19-0.62).

Mean PAE: WT 6.40 -> variant 6.33 (delta -0.07). Local PAE at residue 28 (row+column mean): WT 4.33 -> variant 4.32.

| Residue | 25 | 26 | 27 | 28 | 29 | 30 | 31 |
|---|---|---|---|---|---|---|---|
| pLDDT WT (mean of 5) | 94.3 | 95.7 | 94.7 | 95.5 | 96.0 | 95.7 | 96.2 |
| pLDDT variant (mean of 5) | 94.5 | 96.0 | 95.0 | 96.1 | 96.2 | 95.9 | 96.4 |

**Verdict: within seed noise.** (Rule: variant mean RMSD must exceed the maximum WT-vs-WT value.)

### E62G (residue 62: E->G, verified in all 5 WT and 5 variant CIFs)

Reference: WT model 0 (top ranking_score).

| Metric | Variant vs WT-best, mean (range) | WT-vs-WT noise, 10 pairs, mean (range) |
|---|---|---|
| Whole-chain CA RMSD (A) | 0.60 (0.54-0.68) | 0.82 (0.51-1.21) |
| CA RMSD residues 57-67 (A) | 0.30 (0.26-0.35) | 0.88 (0.24-1.60) |

WT-vs-WT noise vs the best WT model only (4 pairs): CA 0.69 (0.51-0.83), local 0.65 (0.24-1.57).

Mean PAE: WT 6.40 -> variant 6.40 (delta +0.00). Local PAE at residue 62 (row+column mean): WT 5.51 -> variant 5.76.

| Residue | 59 | 60 | 61 | 62 | 63 | 64 | 65 |
|---|---|---|---|---|---|---|---|
| pLDDT WT (mean of 5) | 97.8 | 97.6 | 97.2 | 94.0 | 82.7 | 86.7 | 93.6 |
| pLDDT variant (mean of 5) | 97.8 | 97.5 | 97.0 | 93.0 | 82.6 | 87.9 | 93.6 |

**Verdict: within seed noise.** (Rule: variant mean RMSD must exceed the maximum WT-vs-WT value.)

## Caveats

- AlphaFold 3 does not predict the stability or binding effect of point mutations; small RMSD/pLDDT shifts are not evidence of functional impact.
- n=5 models per condition, single seed each; the noise floor is a rough estimate.
- Both variants are missense changes in a mostly disordered/loop-rich small protein; interpret with wet-lab data.
