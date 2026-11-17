# Novelty check (v1 bar): top-10 loci, deep notes

Checked 2026-11-17. Research use only; AlphaGenome scores are predictions, not validation. Full 39-lead table: `novelty_all_loci.csv`; verified references: `references.bib`.

## Tier rule as applied (decision needed from the owner)
- Tier 1 = the variant itself is a genome-wide-significant keloid hit in the GWAS Catalog, or a published fine-mapping/functional paper names it.
- Tier 2 = the locus is a published keloid hit, but this variant (an LD proxy of the lead) and its predicted mechanism are not named.
- Tier 3 = unreported variant and gene, no skin eQTL. **None assigned**: every one of the 39 leads is a published Catalog keloid hit, so no locus is unreported.
- Applying the brief literally (an LD proxy at r2>=0.8 of a GWS hit makes a variant "already known") would make nearly every variant tier 1, since every candidate was selected at r2>=0.8 with its lead. I therefore report the strict variant-level tier in `tier` and the literal result in `tier_ld_rule_literal`. Both are in the CSV.
- Light-pass rows (29 leads) used GWAS Catalog + GTEx only; their tier is provisional (`tier_confidence`).

## Tier counts (all 39 leads)
Tier 1: 15, Tier 2: 24, Tier 3: 0 (4 of the tier-1 leads were not scored by AlphaGenome: rs527569697, rs549023067, rs686722, rs76024540). Literal-LD-rule tier 1: 39.

## Top-10 summary
| lead | variant | tier | literal tier | GTEx skin eQTL | min p |
|---|---|---|---|---|---|
| rs6726716 | rs7579101 | 2 | 1 | y | 4.39e-05 |
| rs6476838 | rs2068979 | 2 | 1 | n |  |
| rs140716753 | rs140716753 | 1 | 1 | n |  |
| rs921721 | rs921721 | 1 | 1 | n |  |
| rs646315 | rs940189 | 2 | 1 | y | 2.15e-06 |
| rs1511412 | rs360710 | 2 | 1 | n |  |
| rs13051336 | rs2026232 | 2 | 1 | n |  |
| rs2832056 | rs2049878 | 2 | 1 | y | 1.62e-10 |
| rs35383942 | rs35383942 | 1 | 1 | y | 3.76e-06 |
| rs12568930 | rs12568930 | 1 | 1 | n |  |

Tier 2 (known locus, new variant/mechanism): rs7579101, rs2068979, rs940189, rs360710, rs2026232, rs2049878. Strongest new-mechanism leads: rs2049878 (skin and fibroblast eQTL for N6AMT1/HEMK2), rs940189 (skin eQTL for FOXL2), rs7579101 (BABAM2 eQTL, weak). Closest to tier 3: rs2026232 (LINC03138, no gene literature, no eQTL). rs360710 has no AlphaGenome skin accessibility.

## Per-variant notes
### rs6726716 -> rs7579101  (rank_score 0.998) | Tier 2 | GTEx skin eQTL: y (min p 4.39e-05)
- **Tier basis:** LD proxy r2_EUR=0.90 / r2_EAS=0.99 of Catalog GWS lead rs6726716; proxy not reported; BABAM2 skin eQTL not reported
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: lead rs6726716 (BABAM2) is a GWAS Catalog keloid hit, p=8e-9, Greene 2025. rs7579101 itself: 0 GWAS Catalog associations, 0 rsID hits in PubMed and in Europe PMC open-access full text; ClinVar 0. BABAM2: no keloid/scar/fibrosis paper found (PubMed 'BABAM2 OR BRE' AND keloid/scar/fibrosis returned only unrelated hits). NEW: GTEx v10 not-sun-exposed skin BABAM2 eQTL for rs7579101 p=4.4e-05 (lead rs6726716 p=2.0e-04; not significant in v8 or sun-exposed skin). AlphaGenome: non-coding transcript exon, ATAC+DNASE+CHIP_H skin-accessible, max|effect| 4.7 (prediction).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); Europe PMC rsID full-text search (0 hits); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs7579101; GTEx https://gtexportal.org/home/snp/rs7579101; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs7579101' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs7579101

### rs6476838 -> rs2068979  (rank_score 0.992) | Tier 2 | GTEx skin eQTL: n
- **Tier basis:** LD proxy r2=1.0 (EUR and EAS) of Catalog GWS lead rs6476838; proxy not reported; no skin eQTL
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: lead rs6476838 (GLIS3) Catalog keloid hit p=1e-10, Greene 2025. rs2068979: 0 Catalog associations, 0 rsID hits in PubMed/Europe PMC, ClinVar 0, no GTEx skin or cultured-fibroblast eQTL (v8, v10). GLIS3 has no keloid/scar paper in PubMed (0 hits); fibrosis context only outside skin: Pokatayev 2026 (Nature) identifies GLIS3 as a transcription factor governing inflammatory/fibrotic gene programs in intestinal fibroblasts (IBD), Glis3 fibroblast deletion alleviated colitis in mice. Open Targets keloid association list ranks GLIS3 22nd (score 0.283). AlphaGenome: intron, only CHIP_H skin-accessible, max|effect| 1.5 (prediction).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); PMID:41501466 (Pokatayev 2026, doi:10.1038/s41586-025-09907-x); Open Targets Platform API get-associated-targets MONDO_0005348; GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs2068979; GTEx https://gtexportal.org/home/snp/rs2068979; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs2068979' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs2068979

### rs140716753 -> rs140716753  (rank_score 0.985) | Tier 1 | GTEx skin eQTL: n
- **Tier basis:** variant itself is a Catalog GWS keloid hit (Verma 2024; Greene 2025) AND named putatively causal in Greene fine-mapping (African ancestry, JAG1)
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** ALREADY KNOWN as association and fine-mapped causal candidate: Catalog keloid p=4e-12 (Verma 2024, VA MVP, PheCode 701.4) and p=2e-14 (Greene 2025 multi-ancestry and African analyses). Greene full text: 'putatively causal SNPs included ... rs76024540 (JAG1), and rs140716753 (JAG1) in Africans' (SuSiEx, PIP>0.9 set). Not in ClinVar (0). No functional paper found. Greene's JAG1 label is a mapping: rs140716753 (chr20:11,125,355) lies ~451 kb from the JAG1 TSS (chr20:10,674,396, minus strand; Ensembl), nearest annotated genes are lncRNAs. JAG1/Notch in keloid is published independently of this SNP (Syed 2012; Yu 2023; Yang 2026). No GTEx skin/fibroblast eQTL (v8, v10). AlphaGenome-predicted regulatory effect (regulatory_region_variant, ATAC+DNASE+CHIP_H, max|effect| 1.69) is not named in any source found (prediction).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); PMID:39024449 (Verma 2024, Catalog accession GCST90480466 only); PMID:22985040 (Syed 2012); PMID:36734859 (Yu 2023); PMID:42769688 (Yang 2026); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs140716753; GTEx https://gtexportal.org/home/snp/rs140716753; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs140716753' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs140716753

### rs921721 -> rs921721  (rank_score 0.971) | Tier 1 | GTEx skin eQTL: n
- **Tier basis:** variant itself is the Catalog GWS keloid lead (Greene 2025, p=2e-8)
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** ALREADY KNOWN as association: Catalog keloid p=2e-8, Greene 2025 (multi-ancestry; mapped TRIB1AL). No fine-mapping/functional paper naming it found (Europe PMC rsID: Greene only; PubMed 0). Gene context is ambiguous: Ensembl TRIB1AL (lncRNA, overlaps variant), dbSNP annotates LINC02964; nearest protein-coding TRIB1 TSS is 92 kb away (chr8:125,430,152). No keloid/scar/fibrosis literature for TRIB1/TRIB1AL found (PubMed 0 for keloid AND TRIB1). No GTEx skin or fibroblast eQTL (v8, v10); ClinVar 0. AlphaGenome: intron, ATAC+DNASE+CHIP_H, max|effect| 1.16 (prediction); mechanism unreported.
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs921721; GTEx https://gtexportal.org/home/snp/rs921721; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs921721' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs921721

### rs646315 -> rs940189  (rank_score 0.946) | Tier 2 | GTEx skin eQTL: y (min p 2.15e-06)
- **Tier basis:** LD proxy (r2_EAS=0.90 to Catalog lead rs646315; r2 0.99 EUR/0.98 EAS to Nakashima 2010 lead rs940187) not itself reported; skin FOXL2 eQTL not reported
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: rs646315 is a Catalog keloid hit (Ishigaki 2020 p=4e-12; Sakaue 2021 p=2e-9/5e-10; BioBank Japan, East Asian) and Nakashima 2010 reported rs940187 (p=1.8e-13) in this 3q22.3 block, replicated in Chinese Han (Lu 2018). rs940189 is r2=0.994 (EUR) / 0.976 (EAS) / 0.71 (AFR) with rs940187 (Ensembl LD REST, 1000G P3), so it tags the known signal but is not named in any source found (rsID search 0 hits; 0 Catalog associations). Lu 2018 reported only trans eQTLs for rs940187 (SLC7A9, LEMD3; p<1e-4 in annotation datasets). NEW: rs940189 is a GTEx skin eQTL for FOXL2: v10 sun-exposed p=8.0e-06, not-sun-exposed p=2.1e-06 (v8: 8.9e-05, 3.5e-05); also FOXL2NB (v10 not-sun p=2.2e-05). FOXL2 is the gene Nakashima/Lu 2015 attach to the separate block rs1511412; rs940189 is only r2=0.42 (EUR) / 0.45 (EAS) with rs1511412, so the two blocks may converge on FOXL2 regulation (untested hypothesis). Catalog-mapped genes BPESC1/MRPS22: no keloid paper (PubMed 0). AlphaGenome: intron, ATAC+DNASE+CHIP_H, max|effect| 0.56 (prediction).
- **Sources:** PMID:20711176 (Nakashima 2010); PMID:29511450 (Lu 2018); PMID:32514122 (Ishigaki 2020, Catalog accession only); PMID:34594039 (Sakaue 2021, Catalog accession only); PMID:25932232 (Lu 2015); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs940189; GTEx https://gtexportal.org/home/snp/rs940189; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs940189' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs940189

### rs1511412 -> rs360710  (rank_score 0.94) | Tier 2 | GTEx skin eQTL: n
- **Tier basis:** LD proxy r2_EUR=0.83 of Catalog GWS lead rs1511412; proxy not reported; AlphaGenome shows no skin accessibility
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: lead rs1511412 Catalog keloid hit p=2e-13 (Nakashima 2010, Japanese; also associated with 'Balding type 1' in Catalog), replicated in Chinese Han meta-analysis (Lu 2015: P<1e-8, 5 Asian populations; Lu 2018 only a trend). rs360710: 0 Catalog associations, 0 rsID hits, ClinVar 0, no GTEx skin/fibroblast eQTL for rs360710 or rs1511412 (v8, v10). FOXL2 is the gene named for rs1511412 ('SNP rs1511412 in FOXL2 gene'); rs360710 lies 27.5 kb from the FOXL2 TSS (chr3:138,947,137, minus strand). Caution: AlphaGenome rank is high but skin_accessible=False, so the regulatory prediction lacks skin chromatin support. Mechanism not named (prediction).
- **Sources:** PMID:20711176 (Nakashima 2010); PMID:25932232 (Lu 2015); PMID:29511450 (Lu 2018); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs360710; GTEx https://gtexportal.org/home/snp/rs360710; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs360710' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs360710

### rs13051336 -> rs2026232  (rank_score 0.937) | Tier 2 | GTEx skin eQTL: n
- **Tier basis:** LD proxy r2_EUR=0.81 of Catalog GWS lead rs13051336; no gene literature, no skin eQTL (closest to tier 3)
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: lead rs13051336 (LINC03138) Catalog keloid hit p=2e-14, Greene 2025. rs2026232: 0 Catalog associations, 0 rsID hits in PubMed/Europe PMC, ClinVar 0, no GTEx skin/fibroblast eQTL (v8, v10). LINC03138: no paper in PubMed (keloid AND LINC03138 = 0); no protein-coding gene within 190 kb. Closest to tier 3 of the top ten but stays tier 2 because the locus is a published hit. AlphaGenome: intron, ATAC+DNASE+CHIP_H, max|effect| 1.06 (prediction).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs2026232; GTEx https://gtexportal.org/home/snp/rs2026232; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs2026232' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs2026232

### rs2832056 -> rs2049878  (rank_score 0.915) | Tier 2 | GTEx skin eQTL: y (min p 1.62e-10)
- **Tier basis:** LD proxy r2_EUR=0.83 of Catalog GWS lead rs2832056; proxy not reported; strong skin/fibroblast eQTL for N6AMT1(HEMK2) not reported in sources read
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** Locus known: lead rs2832056 (LINC03138) Catalog keloid hit p=1e-9, Greene 2025. rs2049878: 0 Catalog associations, 0 rsID hits, ClinVar 0. NEW: GTEx v10 skin eQTL for N6AMT1 (=HEMK2): sun-exposed p=1.6e-10, cultured fibroblasts p=2.1e-9 (v8 sun p=3.6e-05); also RPL23P2, CCT8, LINC00161. Lead rs2832056 is also a skin eQTL for N6AMT1 (v10 sun p=4.7e-10; fibroblasts p=1.4e-11). Open Targets lists HEMK2 as the 3rd-ranked keloid-associated target (score 0.446; evidence type not examined), but no keloid/scar/fibrosis paper on HEMK2/N6AMT1 was found (PubMed 0) and Greene's main text does not discuss it (supplementary tables not read). HEMK2 TSS is ~197 kb from the variant (chr21:28,885,384, minus strand). AlphaGenome: intron, ATAC+DNASE+CHIP_H, max|effect| 0.68 (prediction).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); Open Targets Platform API get-associated-targets MONDO_0005348 (HEMK2 rank 3); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs2049878; GTEx https://gtexportal.org/home/snp/rs2049878; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs2049878' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs2049878

### rs35383942 -> rs35383942  (rank_score 0.915) | Tier 1 | GTEx skin eQTL: y (min p 3.76e-06)
- **Tier basis:** variant itself is a Catalog GWS keloid hit (Greene 2025; Dand 2025 'excessive scarring') AND named putatively causal in Greene fine-mapping
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** ALREADY KNOWN: Catalog keloid p=2e-25 (Greene 2025) and excessive scarring p=7e-11 (Dand 2025, Catalog record only). Greene full text names rs35383942 (PHLDA3) as European putatively causal SNP, T allele OR 1.49, p=1.1e-23, highest CADD (24.3), RegulomeDB 2a, active TSS in dermal fibroblasts (E126); S-PrediXcan: higher predicted PHLDA3 in not-sun-exposed skin associated with keloid risk (p=1.7e-14). Pena 2026 review states it is a credible variant for hypertrophic scar in FinnGen (secondhand). GTEx confirms rs35383942 as skin eQTL for PHLDA3: v10 sun p=4.8e-06, not-sun p=1.4e-05 (v8: 9.1e-05, 3.8e-06). Missense variant; many Catalog pQTL associations (TMPRSS5, ERBB4, GFRA3, MIA). No functional wet-lab paper found. AlphaGenome: missense, ATAC only, max|effect| 0.32 (prediction; weak).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); PMID:39746571 (Dand 2025, Catalog accession GCST90483368 only); PMID:41320625 (Pena 2026 review, secondhand); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs35383942; GTEx https://gtexportal.org/home/snp/rs35383942; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs35383942' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs35383942

### rs12568930 -> rs12568930  (rank_score 0.912) | Tier 1 | GTEx skin eQTL: n
- **Tier basis:** variant itself is a Catalog GWS keloid hit (Greene 2025, p=6e-9)
- **Literal LD-rule tier (r2>=0.8 proxy counts as known):** 1
- **Evidence:** ALREADY KNOWN as association: Catalog keloid p=6e-9, Greene 2025 (mapped PPIAP34;ZBTB40). Also Catalog hits for IBD, bone mineral density, hematocrit/hemoglobin/red-cell traits (not skin). No fine-mapping/functional keloid paper naming it found (Europe PMC rsID: Greene plus unrelated IBD/endometriosis papers). Nearest protein-coding genes: ZBTB40 (TSS 53 kb, chr1:22,428,838) and WNT4 (TSS 232 kb, chr1:22,143,969, minus strand). No keloid paper on ZBTB40 (PubMed 0). WNT4 has independent skin-fibrosis literature: reduced Wnt4 in hypertrophic scar and Wnt4 suppresses TGF-b1 fibroblast-to-myofibroblast transition (Liu 2020); WNT4 downregulation gives Dupuytren-like fibroblast transcriptome (O'Brien 2024); no source links this variant to WNT4. No GTEx skin/fibroblast eQTL (v8, v10). AlphaGenome: regulatory_region_variant, only CHIP_H skin-accessible, max|effect| 0.32 (prediction; weak).
- **Sources:** PMID:40835838 (Greene 2025, doi:10.1038/s41467-025-62945-x); PMID:31776823 (Liu 2020); PMID:39381224 (O'Brien 2024); GWAS Catalog https://www.ebi.ac.uk/gwas/variants/rs12568930; GTEx https://gtexportal.org/home/snp/rs12568930; Ensembl REST (VEP/LD) https://rest.ensembl.org; ClinVar esearch 'rs12568930' = 0 records; dbSNP https://www.ncbi.nlm.nih.gov/snp/rs12568930

## Source conflicts and caveats
- Gene annotation differs by source: rs921721 is TRIB1AL (Ensembl) vs LINC02964 (dbSNP); rs2049878 is LINC03138 (Ensembl VEP/Catalog mapping) vs N6AMT1 (dbSNP) with HEMK2 as the eQTL target.
- GTEx releases differ: BABAM2 eQTL for rs7579101 and rs6726716 appears only in v10 (not v8); the brief asked for v8, I ran v8 and v10 (the project GTEx skill defaults to v10).
- LD is population dependent: rs940189 has no r2_EUR in the ranking file; rs646315 r2 is 0.90 (EAS) but 0.55 (EUR) in Ensembl 1000G; it is 0.99 (EUR) with rs940187.
- Greene 2025 labels rs140716753 "JAG1", but it is ~451 kb from the JAG1 TSS; this is a mapping label, not demonstrated regulation.
- The brief says 35 ranked loci; the file has 39 leads (35 scored, 4 not scored). All 39 are in the CSV.

## Limitations
- rsID literature search is weak: PubMed does not index rsIDs, and the Europe PMC skill is restricted to open-access text. Paywalled functional papers could be missed.
- Greene 2025 supplementary tables (fine-mapping PIPs, S-PrediXcan gene lists, FUMA gene mapping) were not read; they may name variants or genes absent from the main text. Dand 2025, Sakaue 2021, Ishigaki 2020 and Verma 2024 were used only as GWAS Catalog records, not read.
- GTEx check is exact-variant only (no LD-proxy or locus-wide eQTL scan), two skin tissues plus cultured fibroblasts, significant eQTLs only. Absence is not evidence of no effect. GTEx is non-disease tissue, not keloid.
- Open Targets returned only the top 25 of 1,500 keloid-associated targets (HEMK2 rank 3, JAG1 rank 7, GLIS3 rank 22); evidence type behind each score was not examined.
- Light-pass tiers have no literature search; all absence-of-literature statements are limited to PubMed/Europe PMC queries listed above.

## Databases and tools used (accessed 2026-11-17)
GWAS Catalog REST v2 (https://www.ebi.ac.uk/gwas/rest/api/v2); GTEx Portal API v2 (https://gtexportal.org/api/v2; gtex_v8 and gtex_v10; used via project gtex_database skill terms, license notice on file); Ensembl REST (https://rest.ensembl.org: variation, VEP, LD, lookup); ClinVar and dbSNP (project skills, NCBI E-utilities); PubMed (project pubmed_database skill); Europe PMC (project literature_search_europepmc skill, open access only; papers used: PMC12368108, PMC7618695); Open Targets Platform (project opentargets_database skill, disease MONDO_0005348).

## References (consistent format: authors. title. journal. year. PMID/DOI)
- [40835838] Greene CA, Hampton G, Jaworski J, Shuey MM, Khan A, Luo Y, et al.. Multi-ancestry meta-analysis of keloids uncovers novel susceptibility loci in diverse populations. Nature communications. 2025. PMID:40835838; doi:10.1038/s41467-025-62945-x. *Verification: full text read (Europe PMC PMC12368108).*
- [20711176] Nakashima M, Chung S, Takahashi A, Kamatani N, Kawaguchi T, Tsunoda T, et al.. A genome-wide association study identifies four susceptibility loci for keloid in the Japanese population. Nature genetics. 2010. PMID:20711176; doi:10.1038/ng.645. *Verification: abstract read (PubMed).*
- [23667473] Zhu F, Wu B, Li P, Wang J, Tang H, Liu Y, et al.. Association study confirmed susceptibility loci with keloid in the Chinese Han population. PloS one. 2013. PMID:23667473; doi:10.1371/journal.pone.0062377. *Verification: abstract read (PubMed).*
- [29511450] Lu MZ, Ang QQ, Zhang X, Zhang LF, Yao XH, Lv H, et al.. Genomic risk variants at 3q22.3 are associated with keloids in a Chinese Han population. American journal of translational research. 2018. PMID:29511450. *Verification: abstract read (PubMed).*
- [25932232] Lu W, Zheng X, Liu S, Ding M, Xie J, Yao X, et al.. SNP rs1511412 in FOXL2 gene as a risk factor for keloid by meta analysis. International journal of clinical and experimental medicine. 2015. PMID:25932232. *Verification: abstract read (PubMed).*
- [36689521] Deng CC, Zhang LX, Xu XY, Zhu DH, Cheng Q, Ma S, et al.. Risk single-nucleotide polymorphism-mediated enhancer-promoter interaction drives keloids through long noncoding RNA down expressed in keloids. The British journal of dermatology. 2023. PMID:36689521; doi:10.1093/bjd/ljac025. *Verification: abstract read (PubMed).*
- [41501466] Pokatayev V, Jaiswal A, Shih AR, Segerstolpe Å, Li B, Creasey EA, et al.. Bidirectional CRISPR screens decode a GLIS3-dependent fibrotic cell circuit. Nature. 2026. PMID:41501466; doi:10.1038/s41586-025-09907-x. *Verification: abstract read (PubMed).*
- [22985040] Syed F, Bayat A. Notch signaling pathway in keloid disease: enhanced fibroblast activity in a Jagged-1 peptide-dependent manner in lesional vs. extralesional fibroblasts. Wound repair and regeneration : official publication of the Wound Healing Society [and] the European Tissue Repair Society. 2012. PMID:22985040; doi:10.1111/j.1524-475X.2012.00823.x. *Verification: abstract read (PubMed).*
- [36734859] Yu Y, Dong Y, Deng B, Yang T. IncRNA MIAT Accelerates Keloid Formation by miR-411-5p/JAG1 Axis. Critical reviews in eukaryotic gene expression. 2023. PMID:36734859; doi:10.1615/CritRevEukaryotGeneExpr.2022044734. *Verification: abstract read (PubMed).*
- [42769688] Yang L, Mo W, Cao W, Mi L, Wu X. Integrative Transcriptomic and Machine Learning Analysis Identifies JAG1 and NR2F2 in Keloid. Clinical, cosmetic and investigational dermatology. 2026. PMID:42769688; doi:10.2147/CCID.S625100. *Verification: abstract read (PubMed).*
- [31776823] Liu J, Zhao B, Zhu H, Pan Q, Cai M, Bai X, et al.. Wnt4 negatively regulates the TGF-β1-induced human dermal fibroblast-to-myofibroblast transition via targeting Smad3 and ERK. Cell and tissue research. 2020. PMID:31776823; doi:10.1007/s00441-019-03110-x. *Verification: abstract read (PubMed).*
- [39381224] O'Brien A, Stevenson A, Barrett L, Lawler NB, Hortin N, Deng Z, et al.. Reduced WNT4 expression in normal skin fibroblasts leads to 'Dupuytren-like' changes in the transcriptome. Heliyon. 2024. PMID:39381224; doi:10.1016/j.heliyon.2024.e38016. *Verification: abstract read (PubMed).*
- [41320625] Peña OA, Cañadas-Garre M, Collins I, Timpson NJ, Martin P. Harnessing population genetics and animal models to uncover the genetics of skin scarring. Journal of plastic, reconstructive & aesthetic surgery : JPRAS. 2026. PMID:41320625; doi:10.1016/j.bjps.2025.11.010. *Verification: passage read in Europe PMC full text PMC7618695; FinnGen claim is secondhand.*
- [39746571] Dand N, Ung CY, Saklatvala JR, Simpson MA, Barker JN, Shaw TJ, et al.. GWAS Meta-Analysis Identifies Susceptibility Loci for Keloids and Hypertrophic Scarring in Europeans. The Journal of investigative dermatology. 2025. PMID:39746571; doi:10.1016/j.jid.2024.12.011. *Verification: existence verified in PubMed; abstract/full text NOT read; cited only for GWAS Catalog accession GCST90483368.*
- [34594039] Sakaue S, Kanai M, Tanigawa Y, Karjalainen J, Kurki M, Koshiba S, et al.. A cross-population atlas of genetic associations for 220 human phenotypes. Nature genetics. 2021. PMID:34594039; doi:10.1038/s41588-021-00931-x. *Verification: existence verified in PubMed; paper not read; cited only for GWAS Catalog accessions GCST90018874/GCST90018654.*
- [32514122] Ishigaki K, Akiyama M, Kanai M, Takahashi A, Kawakami E, Sugishita H, et al.. Large-scale genome-wide association study in a Japanese population identifies novel susceptibility loci across different diseases. Nature genetics. 2020. PMID:32514122; doi:10.1038/s41588-020-0640-3. *Verification: existence verified in PubMed; paper not read; cited only for GWAS Catalog accession GCST90013715.*
- [39024449] Verma A, Huffman JE, Rodriguez A, Conery M, Liu M, Ho YL, et al.. Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program. Science (New York, N.Y.). 2024. PMID:39024449; doi:10.1126/science.adj1182. *Verification: existence verified in PubMed; paper not read; cited only for GWAS Catalog accessions GCST90480466/GCST90476199.*
