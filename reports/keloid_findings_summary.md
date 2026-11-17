# Keloid computational findings (research only, not clinical)
Generated 2026-11-16. AlphaGenome 2 (API) + AlphaFold 3 (AlphaFold Server). Predictions, not validation.

## 1. Regulatory scoring (AlphaGenome 2, 6 skin cell conditions)
50 GWAS/literature variants, 70 variant-allele pairs. Rank score = average of per-track percentile ranks; it is an
ad hoc summary, **not** a validated pathogenicity score. Top 10:

| rs_id       | alt   | gene              | consequence               |   rank_score |   p_value | in_gwas_catalog_keloid   |
|:------------|:------|:------------------|:--------------------------|-------------:|----------:|:-------------------------|
| rs76024540  | T     | SLC67A1-AS        | intron_variant            |        1     |     9e-16 | True                     |
| rs140716753 | A     | LINC02871;RPS11P1 | regulatory_region_variant |        0.973 |     4e-12 | True                     |
| rs921721    | T     | TRIB1AL           | intron_variant            |        0.967 |     2e-08 | True                     |
| rs35383942  | T     | PHLDA3            | missense_variant          |        0.927 |     2e-25 | True                     |
| rs12568930  | C     | PPIAP34;ZBTB40    | regulatory_region_variant |        0.914 |     6e-09 | True                     |
| rs75826502  | C     | BPESC1;MRPS22     | intron_variant            |        0.879 |     1e-09 | True                     |
| rs4239705   | G     | LINC02871;RPS11P1 | intron_variant            |        0.864 |     5e-15 | True                     |
| rs1442440   | C     | nan               | intergenic_variant        |        0.847 |   nan     | False                    |
| rs34647667  | G     | CORO2B;ITGA11     | intron_variant            |        0.806 |     1e-41 | True                     |
| rs747722    | T     | nan               | intron_variant            |        0.795 |   nan     | False                    |

## 2. Structure (AlphaFold 3, wild-type vs variant)
A change counts only if it exceeds the wild-type seed-to-seed noise.

| variant   |   ca_rmsd_mean |   noise_ca_mean |   local_rmsd_mean |   noise_local_mean | verdict           |
|:----------|---------------:|----------------:|------------------:|-------------------:|:------------------|
| R28Q      |          0.677 |           0.816 |             0.344 |              0.42  | within seed noise |
| E62G      |          0.605 |           0.816 |             0.302 |              0.877 | within seed noise |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |
| nan       |        nan     |         nan     |           nan     |            nan     | nan               |

## 3. Agreement with prior evidence
_PENDING: benchmark (reports/benchmark_vs_prior.md not found)_

## 4. Biological context
# Biological context notes: PHLDA3 R28Q / E62G, MYO7A V679I, top-5 AlphaGenome variants

Research-only. Outputs are database look-ups and predictions, not validation. All data retrieved **2026-11-16** (every row below; the date is the retrieval date unless a paper date is given). Raw JSON is in `data/processed/context/`.
Sources used: gnomAD (gnomad_r4, GRCh38), ClinVar (NCBI), InterPro, UniProt, Ensembl VEP (human, GRCh38), GTEx v10, Human Protein Atlas, UCSC API (ENCODE cCRE v3 track), Europe PMC, PubMed. Per the skills, users should check the terms of gnomAD, ClinVar, InterPro, GTEx, HPA, ENCODE (notice files in `.licenses/`).

## 1. Variant identity and AlphaGenome rank

| Variant | GRCh38 | Protein change | AlphaGenome rank (of 70) / score | In GWAS Catalog keloid set? | Source |
|---|---|---|---|---|---|
| rs35383942 | 1:201468704 C>T | PHLDA3 R28Q (Q9Y5J5) | 4 / 0.927 | Yes (Greene 2025 p=2e-25; Dand 2025 p=7e-11) | variant_ranking.csv; keloid_gwas_variants.csv; VEP 2026-11-16 |
| rs192314256 | 1:201468602 T>C | PHLDA3 E62G | 24 / 0.594 | Yes (Sakaue 2021, Ishigaki 2020; p=9e-19 to 3e-20) | same |
| rs35641839 | 11:77174855 G>A | MYO7A V679I (Q13402; position 679 on canonical 2215-aa isoform) | 60 / 0.270 | **No** (seeded from an admixture study) | same |

## 2. gnomAD v4 allele frequency (GRCh38, gnomad_r4; joint = exome+genome)

| Variant | Joint AC/AN | Joint AF | Homozygotes | FAF95 popmax | Source |
|---|---|---|---|---|---|
| rs35383942 (1-201468704-C-T) | 86518 / 1608206 | 5.4% | 2855 | 6.07% (nfe) | gnomAD API 2026-11-16 |
| rs192314256 (1-201468602-T-C) | 662 / 1613812 | 0.041% | 4 | 1.24% (eas) | gnomAD API 2026-11-16 |
| rs35641839 (11-77174855-G-A) | 5194 / 1613706 | 0.32% | 124 | 5.88% (afr) | gnomAD API 2026-11-16 |

Ancestry-specific AF (exome / genome; gnomAD ancestry labels, XX/XY and HGDP/1KG subsets omitted, in raw JSON):

| Ancestry | R28Q (T) | E62G (C) | V679I (A) |
|---|---|---|---|
| nfe (non-Finnish European) | 6.13% / 5.82% | 0.0001% / 0 | 0.014% / 0.016% |
| fin | 11.6% / 12.9% | 0 / 0 | 0.002% / 0 |
| afr | 0.86% / 0.97% | 0 / 0 | **6.06% / 6.00%** |
| eas | 0.010% / 0 | **1.34% / 1.31%** | 0.0025% / 0 |
| sas | 0.70% / 0.73% | 0.005% / 0.041% | 0.015% / 0.021% |
| amr | 2.7% / 3.8% | 0.0045% / 0.013% | 0.35% / 0.63% |
| asj | 3.1% / 3.0% | 0 / 0 | 0 / 0 |
| mid | 2.2% / 2.0% | 0 / 0 | 0.31% / 0.34% |

Gene constraint (gnomAD, 2026-11-16): PHLDA3 pLI 3e-6, LoF o/e 1.58 (upper CI 2.55), missense o/e 1.01 (0.89-1.14), so neither LoF- nor missense-constrained. MYO7A pLI 1e-55, LoF o/e 0.80 (0.71-0.89), missense o/e 0.95 (0.92-0.98); MYO7A disease is mostly recessive, so low pLI is expected.

## 3. ClinVar

| Variant | ClinVar status | Source |
|---|---|---|
| rs35383942 R28Q | No ClinVar record found (rsID query `35383942[RS]` returned 0; gene-wide PHLDA3 list of 43 records has no p.Arg28Gln) | ClinVar 2026-11-16 |
| rs192314256 E62G | VariationID 4343623, NM_012396.5:c.185A>G (p.Glu62Gly): clinical significance "Unknown", review "Unknown" (no clinical assertion) | ClinVar 2026-11-16 |
| rs35641839 V679I | VariationID 43170, NM_000260.4:c.2035G>A: **Benign/Likely benign**, criteria provided, multiple submitters, no conflicts, last evaluated 2026/02/02 (conditions listed: Usher 1, DFNB2, DFNA11) | ClinVar 2026-11-16 |

Other PHLDA3 ClinVar missense records are all VUS (single submitter); none are keloid-related.

## 4. Domain location

| Variant | Residue | UniProt (entry v177 / v255, annotated 2026-11-04) | InterPro (2026-11-16) |
|---|---|---|---|
| R28Q | PHLDA3 Arg28 | Domain "PH" 8-127 (covers nearly the whole 127-aa protein) | Inside PH domain: IPR001849 (9-105), Pfam PF00169 (9-88), SMART SM00233, CDD cd00821 (11-88), PHLDA family IPR042832 (3-122). UniProt lists R28 as natural variant VAR_050515 (rs35383942) |
| E62G | PHLDA3 Glu62 | Same PH domain 8-127 | Same entries (all contain residue 62) |
| V679I | MYO7A Val679 | Domain "Myosin motor" 65-741 | Inside myosin motor domain: IPR036106 / CDD cd01381 (79-729), Pfam PF00063 (67-729), PROSITE PS51456 (65-741). CDD residue annotation places V679 in the **"converter subdomain"** site (residues 672-694+). Not in IQ, MyTH4, FERM or SH3 domains |

Not determined here: whether R28 or E62 contact the phosphoinositide-binding pocket or an AKT interface (needs the structure step).

## 5. Ensembl VEP (GRCh38; SIFT/PolyPhen), re-verified

| Variant | SIFT | PolyPhen | Note | Source |
|---|---|---|---|---|
| R28Q | deleterious (0.01), all 7 transcripts | benign (0.32) | Matches prior notes | Ensembl REST VEP 2026-11-16 |
| E62G | deleterious (0) | possibly_damaging (0.678) | Matches prior notes | same |
| V679I | **Mixed across transcripts**: tolerated_low_confidence (0.05-0.06) on ENST00000409709, 00000660626, 00000669443, 00000962571, 00000973085, 00000973090, 00001141039; deleterious_low_confidence (0-0.03) on ENST00000409619, 00000458637, 00000973086-89 | benign (0-0.17) where scored | Prior note "tolerated" is only partly right: all calls are low-confidence and sit at the 0.05 cutoff | same |

## 6. Expression

| Gene | GTEx v10 median TPM (skin, fibroblast) | GTEx rank among 54 tissues | HPA protein (IHC) | Source |
|---|---|---|---|---|
| PHLDA3 | Cultured fibroblasts 87.6; skin sun-exposed 69.9; skin not sun-exposed 61.4 | 9th, 20th, 21st (top: tibial nerve 283.8) | Skin 1 medium, Skin 2 medium; soft tissue low; adipose not detected; RNA "low tissue specificity" | GTEx API; HPA API 2026-11-16 |
| MYO7A | Cultured fibroblasts **0.042**; skin sun-exposed 0.94; skin not sun-exposed 0.89 | 53rd, 46th, 49th (top: testis 43.6) | Skin 1 medium, Skin 2 high; soft tissue not detected; smooth muscle not detected | same |

HPA flags RNA-protein agreement as low for both genes ("approved" reliability for PHLDA3, "uncertain" for MYO7A), and the HPA wrapper does not give cell-type resolution inside skin, so fibroblast protein levels are not established. MYO7A RNA is essentially absent in cultured fibroblasts, which weakens a fibroblast-intrinsic mechanism.

## 7. GTEx eQTL (skin sun-exposed, skin not sun-exposed, cultured fibroblasts; GENCODE v10 gene IDs)

| Variant | Result | Source |
|---|---|---|
| rs35383942 | Significant eQTL for **PHLDA3**: sun-exposed skin p=4.8e-6, NES +0.283; non-sun-exposed skin p=1.4e-5, NES +0.263 (T allele raises expression). Not significant in fibroblasts. Stronger PHLDA3 skin eQTLs exist (best p=1.5e-26, rs4915527), so this variant is not the lead eQTL and may only tag one | GTEx 2026-11-16 |
| rs192314256 | Not among PHLDA3 significant eQTLs in the 3 tissues (rare in GTEx donors) | same |
| rs35641839 | Not among MYO7A significant eQTLs in the 3 tissues | same |

## 8. Known biology and literature

| Item | Finding | Source (date) |
|---|---|---|
| PHLDA3 function | p53 target gene encoding a PH-domain-only protein; competes with Akt PH domain for membrane lipids, blocking Akt membrane translocation/activation; loss raises Akt activity and lowers p53-dependent apoptosis; locus often lost in primary lung cancer | Kawase et al., Cell 2009 (PMID 19203586, PubMed abstract); UniProt Q9Y5J5 |
| PHLDA3 in cancer | Reviewed as p53-regulated tumor suppressor linked to PI3K/AKT, EMT, Wnt/beta-catenin | Kamel, Krishnaraj, Ohki, Cancers 2025 (PMID 40227573, abstract) |
| PHLDA3 and keloid | Keloid GWAS lead locus: rs192314256 in East Asians (BBJ), rs35383942 in Europeans (OR 1.49, p=1.1e-23); rs35383942 CADD 24.3 (max of the lead SNPs), RegulomeDB 2a, in an active TSS state in dermal fibroblast epigenome E126. TWAS: higher predicted PHLDA3 expression in non-sun-exposed skin associated with higher keloid risk (OR 2.39, p=1.7e-14). Not significant in the African ancestry analysis (p=0.005, MAF <1%) | Greene et al., Nat Commun 2025 (PMID 40835838, full text) |
| PHLDA3 and scarring (other) | PHLDA3-rs35383942 listed as credible variant for hypertrophic scar in FinnGen; the review notes these genes were not replicated in UK Biobank | Pena et al., Br J Plast Surg 2026 (PMID 41320625, full text) |
| PHLDA3 functional keloid/fibrosis papers | None found in PubMed or Europe PMC (abstracts): "PHLDA3 AND keloid" returned 0 PubMed hits | PubMed, Europe PMC 2026-11-16 |
| MYO7A function | Unconventional myosin (actin motor); retina photoreceptor/RPE, inner ear hair cell stereocilia. Disease: Usher 1B, DFNB2, DFNA11. Not expressed in brain; also kidney, liver, testis, lymphocytes | UniProt Q13402 (2026-11-04 annotation) |
| MYO7A and keloid | Single source: admixture mapping/exome study of 478 African Americans (122 cases, 356 controls); rs35641839 OR 4.71 (2.38-9.32), **p=8.3e-6** (not genome-wide significant), same study also found MYO1E | Velez Edwards et al., Hum Genet 2014 (PMID 25280642, abstract). Cited second-hand by three reviews (PMC11378114, PMC7264387, PMC7618695) |
| MYO7A and fibrosis | No primary paper linking MYO7A to fibroblast function or fibrosis found; keyword hits were reviews that cite Velez Edwards 2014 | Europe PMC 2026-11-16 |

## 9. Top 5 AlphaGenome-ranked variants

| Rank | Variant, alt | AlphaGenome gene / consequence | Nearest genes (Ensembl REST overlap, +/-150 kb) | ENCODE cCRE overlap (UCSC encodeCcreCombined) | GTEx skin/fibroblast eQTL | gnomAD joint AF (ancestry) | GWAS keloid p |
|---|---|---|---|---|---|---|---|
| 1 | rs76024540 T (11:2898878) | SLC67A1-AS, intron | SLC67A1-AS (overlaps; SLC67A1 is GTEx name SLC22A18), SLC67A1 840 bp, CDKN1C 13.1 kb, PHLDA2 29.4 kb, NAP1L4 45.6 kb | EH38E1515835 pELS, CTCF-bound | **PHLDA2 in cultured fibroblasts**: p=3.9e-5, NES -0.273. None for SLC22A18 or CDKN1C | 3.4% (afr 11.9%, nfe 0.04%, eas 0) | 3e-18 (Greene); 9e-16 (Verma) |
| 2 | rs140716753 A (20:11125355) | LINC02871;RPS11P1, regulatory | No protein-coding gene within 150 kb; LINC02871 95.9 kb | EH38E2097882 dELS; Ensembl enhancer ENSR20_8BSG9 | None found (no nearby gene tested, see limits) | 1.0% (afr 3.6%, nfe 0.02%) | 2e-14 (Greene); 4e-12 (Verma) |
| 3 | rs921721 T (8:125522294) | TRIB1AL, intron | TRIB1AL (overlaps), TRIB1 83.3 kb | EH38E2664615 pELS; Ensembl enhancer ENSR8_93SS75 | None in TRIB1 (40 eQTLs retrieved; variant absent) | 59.5% (common in all ancestries) | 2e-8 (Greene) |
| 4 | rs35383942 T (1:201468704) | PHLDA3, missense | PHLDA3 (coding) | EH38E1411712 pELS, CTCF-bound (201468615-201468772) | PHLDA3 skin eQTL, see section 7 | 5.4% (nfe 6.1%, fin 11.6%, afr 0.9%, eas ~0) | 2e-25 (Greene) |
| 5 | rs12568930 C (1:22375738) | PPIAP34;ZBTB40, regulatory | PPIAP34 52 kb, ZBTB40 53.1 kb | EH38E1326968 dELS; Ensembl enhancer ENSR1_9BTRG | None in ZBTB40 (108 eQTLs; variant absent) | 20.7% (afr 32%, fin 7.8%) | 6e-9 (Greene) |

Sources for the table: variant_ranking.csv, keloid_gwas_variants.csv, Ensembl VEP/REST, UCSC API, GTEx API, gnomAD API, all 2026-11-16.

Notes:
- rs76024540 sits in the 11p15.5 gene cluster that includes CDKN1C and PHLDA2 (PHLDA2 is a PHLDA3-family name; the family relationship was not checked in a database here). Greene 2025 text labels rs76024540 and rs140716753 "(JAG1)"; Ensembl and the GWAS Catalog place rs76024540 on chr11 and rs140716753 about 450 kb from JAG1, so the paper's gene label does not match these coordinates (unresolved).
- The ranking `alt` allele differs from the GWAS effect allele for rs12568930 (alt C, GWAS effect T) and rs921721 (alt T, effect A). Ranking uses absolute effect magnitudes, so this affects interpretation of direction only.
- Four of the top five are non-coding (lncRNA intron or enhancer-like cCRE); only rs35383942 is coding. rs76024540 and rs140716753 are African-ancestry-enriched (near zero in nfe/eas in gnomAD), so wet-lab cell sources should match ancestry.

## 10. What weakens the case

- **MYO7A V679I is the weakest candidate.** ClinVar Benign/Likely benign (multiple submitters, no conflicts); 6% allele frequency in African ancestry with 124 gnomAD homozygotes; not in GWAS Catalog; the only evidence is a 122-case admixture study with p=8.3e-6; AlphaGenome rank 60/70; MYO7A RNA is nearly absent in cultured fibroblasts (0.04 TPM); SIFT/PolyPhen are tolerated/benign or low-confidence.
- **PHLDA3 R28Q is common** (5.4% overall, 11.6% Finnish, 2855 homozygotes), so a strongly damaging protein effect is unlikely; PolyPhen is benign and the gene is not missense-constrained (o/e 1.01). The GWAS signal is more plausibly regulatory: the variant sits in a pELS cCRE, is a skin PHLDA3 eQTL (T allele raises expression), and TWAS links higher expression to higher risk. A protein-structure study of R28Q should be framed as exploratory.
- **PHLDA3 E62G** is East Asian-specific (EAS AF 1.3%, essentially absent elsewhere); SIFT 0 and PolyPhen 0.68 favour a protein effect, but AlphaGenome rank is only 24 and ClinVar has no assertion. **Effect direction conflict:** GWAS Catalog gives effect allele C a positive beta (1.83, Sakaue), while Greene 2025 reports "rs192314256-C OR 0.16" (exp(-1.83) = 0.16, so likely an orientation difference). Not resolved here; check the original summary statistics before stating risk vs protective.
- PHLDA3 is LoF-tolerant in gnomAD (o/e 1.58), so loss of function is not obviously harmful at population level.
- eQTL for rs35383942 is not the lead PHLDA3 skin signal, so causality is not established; a nearby variant in LD may drive it (LD not checked).
- Keloid replication: Pena 2026 notes PHLDA3 was not replicated in UK Biobank hypertrophic scar analysis.

## 11. Not found / could not verify

- ENCODE **SCREEN** coordinate search (`screen_api.py search`) returned empty for every query, including the skill's own HBB example and a known accession, so the SCREEN API search appears broken. cCRE overlap was taken instead from the UCSC `encodeCcreCombined` track (ENCODE registry v3, hg38), which is the same registry; cell-type-specific cCRE signal (fibroblast, keratinocyte) was not retrieved.
- GTEx `get-eqtls-in-region` returned empty even for the known rs35383942 positive control, so it was not used. eQTL checks for variants 1, 2, 3, 5 used gene-level queries on nearby protein-coding genes only (SLC22A18, CDKN1C [0 eQTLs returned], PHLDA2, TRIB1, ZBTB40). Not tested: JAG1 or any gene for rs140716753 (no gene within 150 kb); lncRNA targets; other GTEx tissues.
- ClinVar record for R28Q: none found (cannot rule out a record under a different representation).
- HPA: no fibroblast-specific or keratinocyte-specific IHC; subcellular location returned null for both genes; single-cell data not queried.
- Velez Edwards 2014 was read at abstract level only (full text not retrieved); the exact risk allele and frequency in cases/controls are not known.
- Greene 2025 supplementary data (effect-allele tables) not retrieved, so the E62G direction conflict is open.
- Not queried: Open Targets, Reactome, STRING, JASPAR/TF binding, conservation, dbSNP record details, LD. PHLDA2 paralog relationship, 11p15.5 imprinting biology, and the structural position of the residues are not sourced here.
- Europe PMC searches are open-access-only (the wrapper appends `OPEN_ACCESS:y`), so closed-access papers may be missing; PubMed searches returned 0 (PHLDA3+keloid), 1 (MYO7A+keloid).
- Ensembl overlap query for nearest genes used the Ensembl REST API directly (the wrapper has no region-overlap command).


## 5. Limitations
- Most keloid GWAS hits are non-coding and lead SNPs tag nearby causal variants (LD).
- No gold-standard causal variant set exists; scores show concordance, not accuracy.
- AlphaFold does not predict stability or binding effects of single substitutions.
- AlphaFold Server outputs are non-commercial; see results/alphafold/terms_of_use.md.
