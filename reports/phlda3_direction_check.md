# PHLDA3 expression vs keloid risk: direction check

Retrieved 2026-11-16 (Europe PMC full-text XML, GWAS Catalog REST). Research use only; predictions are not validation.

## Verdict: CONFLICT CONFIRMED for rs35383942 (not explained by sign convention); E62G cannot be determined

## Which paper says what
- The TWAS-like result is in **Greene et al., Nat Commun 16:7770 (2025), PMC12368108, doi 10.1038/s41467-025-62945-x** (published 2025-08-20). It uses S-PrediXcan v0.7.1 (GTEx v8 models, 49 tissues) and calls it "genetically predicted gene expression (GPGE)". No z-score is reported for PHLDA3 in the text I read.
- **Dand et al. 2025** (GWAS Catalog GCST90483368, PMID 39746571; Greene ref. 19: J Invest Dermatol 145:1538-1540) is a European-only GWAS (1,779 cases / 434,421 controls; Finland, UK). From the sources I read it has no TWAS or PHLDA3 expression result, so the "Dand TWAS" attribution in our notes is wrong. The expression result belongs to Greene.

## Quotes (Greene 2025, via Europe PMC; text as returned by a summarising fetch tool, so re-check wording against the PDF)
- TWAS direction: "Increased expression of both PHLDA3 in not sun-exposed skin (p = 1.74 x 10-14, OR = 2.39 [1.91-2.99]) and NEDD4 in fibroblasts ... were associated with increased risk of keloids."
- Sign usage elsewhere in the same paper: "increased predicted expression of LINC01705 in fibroblasts was associated with decreased risk of keloids", so direction is stated in plain "predicted expression" terms. The Methods do not define a sign convention.
- Lead SNP (Table 2): rs35383942, effect/risk allele T, ref C, OR 1.50 (1.39-1.61), freq 0.08, p 1.91e-25 ("aligned to reflect increasing risk"). Europeans: "rs35383942-T, p = 1.11e-23, OR = 1.49 [1.38-1.61]". African: "no peak ... (rs35383942-T, p = 0.005)".
- East Asian: "rs192314256-C, p = 8.74e-19, OR = 0.16 [0.11-0.24]" (inconsistent with Sakaue beta +1.83; see e62g_direction_check.md).
- PHLDA3 GPGE was significant in the multi-ancestry and European analyses, not in East Asian or African.
- Cohort ancestry: European, East Asian, African (BioVU, eMERGE, VUMC, MVP, UK Biobank, FinnGen, BBJ); replication in All of Us.
- Not stated in text: the eQTL effect allele and a z-score.

## Independent eQTL check (our own GTEx pull, context_notes.md section 7)
rs35383942 T allele: sun-exposed skin NES +0.283, non-sun-exposed skin NES +0.263 (T raises PHLDA3). T is the risk allele. Greene's direction (higher expression = higher risk) therefore agrees with GTEx (risk allele raises expression). Caveat: rs35383942 is not the lead PHLDA3 skin eQTL (rs4915527, p 1.5e-26); LD was not checked.

## Convention check
A sign-convention mix-up would need Greene's expression-risk sign to be inverted relative to GTEx. It is not: the risk allele T has positive GTEx NES and the TWAS OR is >1 (OR 2.39 for increased expression). Both independent lines give risk allele = higher PHLDA3. So the conflict is with our AlphaGenome prediction (risk allele lowers expression), not a convention artefact. Possible real causes, untested: AlphaGenome tissue/track mismatch (it was scored in fibroblast/keratinocyte, whereas the signal is non-sun-exposed skin; GTEx shows no fibroblast eQTL), AlphaGenome gives a wrong sign for this variant, or the missense/pELS effect acts differently from the net eQTL. I did not re-open the AlphaGenome output to confirm its sign.

## Ishigaki 2020 / Sakaue 2021
Europe PMC lists Ishigaki (PMC7968075), but the full-text fetch failed (HTTP 500), so I found nothing on PHLDA3 there. Sakaue was not re-read; earlier work (e62g_direction_check.md) shows E62G alt C with beta +1.83 (Sakaue) and OR 9.56 (Ishigaki Catalog), with no expression data. E62G expression direction: cannot be determined (rare in GTEx, not an eQTL).

## Not verified
- Exact Greene wording was seen through a summarising tool; the last ~2 KB of the XML and the Dand paper text itself were not read (PubMed blocked).
- Greene S-PrediXcan z-scores/Supplementary tables, and the eQTL effect allele used by the model.
- Whether the AlphaGenome sign in our run is robust.
