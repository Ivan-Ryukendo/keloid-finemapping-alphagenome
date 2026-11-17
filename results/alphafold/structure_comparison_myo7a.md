# MYO7A AlphaFold 3 structure comparison: WT vs V679I (rs35641839)

AlphaFold Server predictions, 5 models per job, 2,215 aa (UniProt Q13402). Research only.
Reference: WT model 0 (top ranking_score). Residue check passed: V679 in all 5 WT models, I679 in all 5 variant models, otherwise identical sequence.

**Domain of residue 679:** Myosin motor domain (UniProt feature 65-741; InterPro IPR036106 79-729, PS51456 65-741, Pfam PF00063 67-729), in the C-terminal part of the motor, before the first IQ motif (745-765). It lies outside the UniProt actin-binding region (632-639). CATH places it in 6.20.240.20 (671-743).

| Metric | Variant vs WT-best, mean (range) | WT-vs-WT noise, 10 pairs, mean (range) | Variant mean > WT-vs-WT max |
|---|---|---|---|
| Whole-chain CA RMSD (A) | 12.51 (7.72-21.29) | 12.44 (8.67-16.84) | no |
| CA RMSD 674-684, global fit (A) | 4.77 (2.31-5.57) | 4.42 (1.32-9.22) | no |
| CA RMSD 654-704, refit on window (A) | 0.23 (0.18-0.32) | 0.25 (0.15-0.40) | no |
| CA RMSD motor domain 65-741, refit on domain (A) | 0.36 (0.27-0.45) | 0.49 (0.29-0.74) | no |

Whole-chain RMSD is dominated by domain/tail movement and is a poor local readout; use the refit rows.

| Confidence | WT (mean of 5) | Variant (mean of 5) |
|---|---|---|
| Mean pLDDT, all residues | 83.4 | 83.5 |
| Fraction of residues pLDDT < 70 | 12.2% | 12.0% |
| Mean pLDDT, residues 654-704 | 85.1 | 85.4 |
| Mean PAE, whole chain (A) | 22.13 | 22.21 |
| Local PAE at residue 679 (row+column mean, A) | 21.16 | 21.14 |

| Residue | 676 | 677 | 678 | 679 | 680 | 681 | 682 |
|---|---|---|---|---|---|---|---|
| pLDDT WT | 86.6 | 87.8 | 91.2 | 90.8 | 87.5 | 90.6 | 90.9 |
| pLDDT variant | 87.3 | 88.3 | 91.5 | 91.3 | 88.0 | 90.8 | 91.4 |

Per model (rank score, mean pLDDT, % residues <70, pLDDT at 679):

| Model | ranking_score | mean pLDDT | % <70 | pLDDT 679 |
|---|---|---|---|---|
| WT 0 | 0.630 | 83.5 | 12.1% | 90.8 |
| WT 1 | 0.620 | 83.4 | 12.1% | 90.8 |
| WT 2 | 0.620 | 83.2 | 13.3% | 90.7 |
| WT 3 | 0.620 | 83.6 | 11.7% | 90.8 |
| WT 4 | 0.620 | 83.5 | 11.9% | 90.7 |
| Variant 0 | 0.620 | 83.3 | 12.7% | 91.1 |
| Variant 1 | 0.620 | 83.6 | 12.1% | 91.3 |
| Variant 2 | 0.620 | 83.8 | 11.7% | 91.4 |
| Variant 3 | 0.610 | 83.2 | 11.6% | 91.5 |
| Variant 4 | 0.610 | 83.6 | 12.0% | 91.2 |

**Verdict: no detectable difference: all variant means are within the WT-vs-WT seed noise.** Rule: variant mean must exceed the maximum WT-vs-WT value for the same metric.

## Caveats

- AlphaFold 3 does not predict the stability or binding effect of point mutations; small RMSD/pLDDT shifts are not evidence of functional impact, and a null result does not show the variant is harmless.
- V679I is a conservative change (Val to Ile) predicted benign/tolerated by VEP.
- n=5 models per condition, single run each; the noise floor is a rough estimate. This is a single-chain prediction without calmodulin/light chains, actin or nucleotide.
