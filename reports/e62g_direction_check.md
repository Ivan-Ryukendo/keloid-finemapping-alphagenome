# PHLDA3 E62G (rs192314256) and R28Q (rs35383942): direction-of-effect check

Retrieved 2026-11-16. Research use only; predictions are not validation.

## Verdict

- **E62G**: the alt allele C (GRCh38 1:201468602 T>C, gAg>gGg, p.Glu62Gly) is associated with INCREASED keloid risk in East Asians. This is shown by two primary datasets that agree: Sakaue 2021 (beta +1.826, p 8.7e-19, harmonised sumstats read directly) and Ishigaki 2020 (OR 9.56, CI 5.91-15.45, Catalog record; I did not read the full paper or its summary statistics).
- The Greene 2025 text "rs192314256-C ... OR = 0.16 (0.11-0.24), p = 8.74e-19" is inconsistent with the above. Its p-value matches the Sakaue/BBJ p-value (8.737e-19) and 1/0.16 = 6.25 = exp(1.826), so it is the reciprocal of the same effect. Most likely the allele label is flipped (0.16 is the OR for T) or it is a typo. I could not verify this, because the Greene East Asian-specific statistics (Supp. Fig. 2) are not deposited in the Catalog and were not retrieved (see Missing).
- The Catalog effect allele C and the paper's stated effect allele C are the same allele by label. The two sources disagree on the direction of that one allele, not on which allele it is.
- **R28Q**: the alt allele T (GRCh38 1:201468704 C>T, cGg>cAg, p.Arg28Gln) is associated with INCREASED risk in Europeans (Greene 2025: OR 1.49, p 1.11e-23; multi-ancestry Table 2 OR 1.50, p 1.91e-25; Catalog GCST90652487). This is consistent across the sources I read. Dand 2025 (GCST90483368) also lists T, p 7e-11, but I did not retrieve its effect size or direction.

## Evidence

| Item | rs192314256 (E62G) | rs35383942 (R28Q) |
|---|---|---|
| GRCh38 position, ref/alt | 1:201468602 T>C | 1:201468704 C>T |
| Ensembl VEP | E/G, gAg/gGg, transcript strand -1 | R/Q, cGg/cAg |
| Allele producing the amino acid change | C (alt) | T (alt) |
| Catalog effect allele | C (Sakaue GCST90018874 and GCST90018654, Ishigaki GCST90013715) | T (Greene GCST90652487, Dand GCST90483368) |
| Catalog effect size | Sakaue beta +1.8259 (CI 1.42-2.23), p 9e-19; Ishigaki OR 9.556 (5.91-15.45), p 3e-20 | Not given for Greene or Dand (beta "-", NR) |
| Sakaue harmonised sumstats (GRCh38) | effect C, other T, beta 1.82592, SE 0.2063, EAF 0.01537, p 8.737e-19 | Variant not present in this file (rare in East Asians) |
| Greene 2025 text | East Asian: "rs192314256-C OR 0.16 (0.11-0.24), p 8.74e-19" | European: "rs35383942-T OR 1.49 (1.38-1.61), p 1.11e-23" |
| Greene Table 2 (multi-ancestry, hg19 coordinates, "aligned to reflect increasing risk") | Not listed | 1:201437832 effect T, ref C, OR 1.50, freq 0.08, p 1.91e-25 |

Frequencies:
- E62G C allele: 1.52% in the study ancestry (BBJ; Sakaue EAF 0.0154, Catalog risk_frequency 0.0152). gnomAD East Asian is about 1.3% (1.34% exome, 1.31% genome), taken from the existing `context_notes.md` gnomAD API pull on 2026-11-16 and not re-queried. These agree, and the Sakaue/Ishigaki effect allele frequency is the C allele, i.e. the rare allele that encodes E62G.
- R28Q T allele: Greene multi-ancestry frequency 0.08. gnomAD overall is 5.4%, non-Finnish European 6.1%, and East Asian about 0 (`context_notes.md`).

Why Sakaue's C is the risk-increasing allele: the harmonised file reports beta for the effect allele (C) against the other allele (T), and the beta is positive. The Catalog states the same value.

## Sources

- GWAS Catalog REST v2: https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id=rs192314256 (Sakaue GCST90018874 and GCST90018654, Ishigaki GCST90013715). Studies: https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90652487 (and /GCST90652488, /GCST90652489, /GCST90018874).
- Sakaue 2021 harmonised summary statistics (GRCh38), read by streaming: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018654/harmonised/GCST90018654.h.tsv.gz. East Asian, N=178,726, PMID 34594039.
- Greene et al., Nat Commun 2025, PMID 40835838, PMC12368108, open access. Full text via https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12368108/fullTextXML. Greene summary statistics are deposited as GCST90652487 (multi-ancestry), GCST90652488 (European) and GCST90652489 (African) only.
- Ensembl REST VEP: https://rest.ensembl.org/vep/human/id/rs192314256 and /rs35383942.
- Ishigaki 2020 (PMID 32514122, PMC7968075) is not open access and I did not read it. Only its Catalog record was used.
- gnomAD values: `reports/context_notes.md` (gnomAD API pull on 2026-11-16).

## Missing

- The Greene East Asian-specific summary statistics (Supplementary Fig. 2 and the study-level BBJ data). Greene's data availability statement points to https://www.pheweb.jp/pheno/Keloid. I could not read the Keloid phenotype row there, and the pheweb.jp variant API returned only a different phenotype for this variant. That is why the 0.16 orientation cannot be confirmed from Greene's own table.
- Ishigaki's paper tables, to confirm the allele orientation of OR 9.56.
- Dand 2025 effect size for rs35383942.
- A Greene Table 2 row for rs192314256 does not exist. The multi-ancestry lead at this locus is rs35383942.

## Wording to use in the report

"The p.Glu62Gly-encoding alt allele (C) at rs192314256 is associated with increased keloid risk in East Asians (BBJ: beta +1.83, OR about 6 to 10, allele frequency 1.5%). The Greene 2025 text reports OR 0.16 for the same allele label with an identical p-value, which is the exact reciprocal and is most likely a labelling error."
