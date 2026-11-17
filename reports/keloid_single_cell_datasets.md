# Public keloid single-cell / spatial datasets

Retrieved 2026-11-16. Sources: NCBI GEO (esearch/esummary on db=gds, GEO acc.cgi text records, ftp.ncbi.nlm.nih.gov/geo/series/*/suppl/ directory listings for file sizes) and Europe PMC (search + fullTextXML). Nothing was downloaded. Research use only; predictions, not validation.

## Datasets

"Normal" = what the GEO record states. Sizes are GEO supplementary file sizes. "Processed" = a ready expression object or matrix is public (beyond per-sample RAW).

| Accession | Paper / year | Tissue, cells | Samples (GEO text) | Processed matrix public? (size) | Notes |
|---|---|---|---|---|---|
| GSE163973 | Deng et al., Nat Commun 2021, PMID 34140509 (author/journal from memory, PMID from GEO) | Dermis; fibroblasts + all cells; keloid vs normal SCAR (not normal skin) | 3 keloid (21,488 cells) + 3 normal scar (19,167 cells) | Yes. Seurat .Rdata: all-cells 1.4G, fibroblast 489M, fib main clusters 445M; cell metadata CSVs 0.3-0.9 MB; RAW.tar 415M | Best fibroblast resource. Seurat objects are over the 200 MB limit; RAW.tar (415M) too. Only the metadata CSV is small. |
| GSE181297 | PMID 35123990 (2022) | Keloid scRNA-seq + Visium spatial | 2 keloid scRNA, 2 keloid ST sections, 2 normal skin ST sections (7 GSM) | Doubtful. `GSE181297_Raw_counts_matrix.txt.gz` is 681 bytes (likely a stub); RAW.tar 316M | Only spatial dataset with normal skin controls. Unverified what is inside RAW.tar. |
| GSE181316 | Schwann cells in keloid, PMIDs 35278628, 36333467 (2022) | Whole skin dissociation | 4 keloid, 3 normal scar, 1 normal skin (8 GSM) | RAW.tar only, 455M | Has normal skin (n=1) and scars. Neurofibroma/healthy comparison data come from external papers. |
| GSE191067 | Midkine immune microenvironment (no PMID on GEO; 2021 submission) | keloid, perilesional, normal scar, normal skin, 3 each | 12 samples, over 100,987 cells | Yes. `all.UMI.matrix.csv.gz` 413M, one matrix | Has true normal skin. Over 200 MB, so not downloaded. |
| GSE270438 | Th17 in keloid, PMID 39872534 (2025) | CD45+ immune cells only | 3 keloid (17,402) + 3 normal scar (23,682) | RAW.tar 239M | Immune only; useless for PHLDA3 in fibroblasts/keratinocytes. |
| GSE220300 | Active vs inactive keloid, PMID 38524141 (2024) | Fibroblasts + endothelial; scRNA + multiplex IF | 11 samples: active/inactive keloid, mature scar, normal scar | RAW.tar 714M | Over limit. |
| GSE335482 | PMID 42433363 (2026) | Keloid infiltrating vs hypercellular zones, chest | 4 patients, paired zones (8 GSM), 128,678 cells | RAW.tar 726M | No normal skin control. |
| GSE282885 | Keloid vs normal tissue, no paper on GEO (2025) | Not stated | 6 samples | RAW.tar 279M | Design text is minimal. Cited by PMC13480444 (J Inflamm Res 2026, keratinocyte EMT). |
| GSE243716 | Hypertrophic scar vs keloid (2024) | Not stated | 2 samples | RAW.tar 234M | Small. |
| GSE274709 / GSE266334 | PIEZO2 fibroblasts, PMID 40742741 (2025) | GSE274709 is bulk RNA-seq (30 samples); GSE266334 is scRNA (4 samples) | see left | 31M / 186M | GSE274709 is NOT single-cell. GSE266334 (186M) is scRNA. |
| GSE307504 | Pathological-scar sc + spatial atlas (2025) | keloid scRNA + spatial, 6 GSM | 6 | RAW.tar 72M; annotated Seurat tar.gz 6.8G | Small RAW. Seurat object far over limit. |
| GSE341250 / GSE325137 | Xenium keloid, Deng et al., Burns Trauma 2026, PMID 42732423 / IL-4R complex, GSE325137 (2026) | Xenium targeted 360-gene panel | Normal skin 1 donor (8,749 cells) + 3 keloid donors (41,264 / 54,831 / 22,489 cells) | RAW.tar not sized | Targeted panel (260-gene skin + 100 custom); whether PHLDA3 is on it is unverified. Normal skin is only one Asian donor. |
| GSE206790 | Acne keloidalis, PMID 38218364 (2024) | Scalp, 1 patient lesional vs non-lesional | 4 GSM | processed bulk CSV 796K only; RAW.tar 157M | Acne keloidalis, not keloid. Raw withheld for privacy. |

Others found in the GEO search but not examined for expression data: GSE237754 (keloid scRNA, 2 samples, mtx 89M, barcodes/genes/tSNE coordinates public), GSE253664 (keloid fibroblasts, xlsx 15M, bulk-like counts), GSE261116 (radiotherapy-treated keloid fibroblasts, 5 samples), GSE175866 (bulk keloid fibroblast, 1.1M), GSE190626 (xlsx 2.1M), GSE303591/303592/293834 (dermal fibroblast subsets, in vitro; keloid relevance unverified).

Literature-only items (no new GEO accession checked): Frontiers in Medicine 2026 multi-dataset scRNA meta-analysis (PMID 42540716, PMC13426371) lists GSE181297, GSE163973, GSE156326, GSE190847, GSE169500, GSE172433, GSE147424 plus HRA000425. Its internal numbers look inconsistent (e.g. text says 8 studies, donor and cell counts differ between text and Table 1), so treat its cohort tables as unverified and use only the accession list. GSE147424 is atopic dermatitis skin per my memory of that accession, not keloid (unverified).

## PHLDA3 in keloid vs normal: what was found

No paper retrieved reports an expression direction (up or down) for PHLDA3 in keloid vs normal skin or fibroblasts. Europe PMC query `PHLDA3 AND keloid` returned 5 records, all genetics or reviews, none with expression data. I read the full text of PMC12368108, PMC7618695, PMC13570645, PMC13426371, PMC13619338, PMC13722046, PMC13683095, PMC13480444; PHLDA3 appeared only in the first two. Verbatim statements, none about expression:

- Greene et al., Nat Commun 2025, PMID 40835838 (PMC12368108): "the other strong result on chromosome 1, consisting of SNPs mapping to PHLDA3, is driven by the European (rs35383942-T, p = 1.11 x 10-23, OR = 1.49 [95% CI 1.38-1.61]) and East Asian analyses (rs192314256-C, p = 8.74 x 10-19, OR = 0.16 [95% CI 0.11-0.24])". Table 2 lists rs35383942 (1:201437832, nearest gene PHLDA3, risk allele freq 0.08, P 1.91 x 10-25). Also: "Six SNPs ... were mapped to PHLDA3, TAB2/SUMO4, LSP1, and NEDD4" (nonsynonymous variants in LD with lead SNPs).
- Review, J Plast Reconstr Aesthet Surg 2026, PMID 41320625 (PMC7618695): "Another gene suggested in a Japanese population from BioBank Japan ... similar Europeans from FinnGen as a marker for keloid formation is PHLDA3 (Pleckstrin Homology Like Domain Family A Member 3)". Same text notes rs35383942 is "a credible variant for hypertrophic scar in the FinnGen study" and that the genes "have not been replicated in the UK BioBank".

Not verified: whether PHLDA3 is differentially expressed in any of the datasets above. This has to be checked in the data itself (e.g. GSE163973 or GSE191067 matrices), not from papers. Other keloid GWAS loci genes (NEDD4, LINC01705 region) were not searched for expression claims.

## Recommendation

1. GSE163973: fibroblast subtypes, keloid vs normal scar, ready Seurat objects, metadata CSVs small. Needs a big download (1.4G or 489M fibroblast-only).
2. GSE191067: single 413M UMI matrix, includes normal skin, perilesional skin, all cell types.
3. GSE181297 (spatial with normal skin ST) or GSE181316 (4 keloid + 1 normal skin + 3 scars, 455M).

## Notes

- Disk space: the C: volume was reported full by the tool during this task. Nothing was written there on purpose. Scratch files are in `scratch_tmp\` (epmc.py, ft.py, ft2.py, geo.py and their .txt outputs), which can be deleted.
- Subagent-free; I did not query Human Protein Atlas or Open Targets for this task.
