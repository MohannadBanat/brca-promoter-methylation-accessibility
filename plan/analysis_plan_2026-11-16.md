# Internal analysis plan committed on 2026-11-16

**Working title:** Promoter DNA Methylation and Chromatin Accessibility as Determinants of Gene Silencing in Basal-Like Breast Cancer: An Integrative Multi-Omic Analysis of TCGA Tumours

**Study type:** secondary, observational, integrative multi-omic analysis of public tumour data (no new samples, no wet lab). Title, abstract and introduction will state this.

**Target journal:** International Journal of Molecular Sciences (MDPI).

**Status at commit:** data downloaded; only file structure and sample counts checked (feasibility gate, Section 3). No association between methylation, accessibility and expression has been computed or viewed.

## 1. Question

In primary breast tumours, how much of a gene's cross-tumour expression variation is explained by promoter DNA methylation, by promoter chromatin accessibility, and by both together? Do genes silenced in basal-like tumours (relative to luminal tumours) differ in which mechanism accounts for the silencing?

Prior work this must be distinguished from: Corces et al. 2018 (TCGA ATAC-seq; peak-to-gene links and ELMER distal-probe analysis); Aine et al. 2025 (TNBC methylation subtypes; ATAC peaks used to filter distal CpGs; 8 cell lines); Sasiain et al. 2026 (CpG cassettes). None reports a per-gene variance decomposition of expression into promoter methylation vs accessibility in tumours with all three assays.

## 2. Data

| Assay | Source | Processing used |
|---|---|---|
| ATAC-seq | Corces et al. 2018, GDC ATACseq-AWG: BRCA-specific peak calls and normalized counts | technical replicates averaged per tumour |
| DNA methylation | GDC harmonized SeSAMe level-3 beta values, one primary-tumour file per case: HM450 (43 tumours) or HM27 (30 tumours) | restricted to probes present on both arrays (HM27 ∩ HM450); hg38 positions and masks from the Zhou lab HM450/HM27 GENCODE v36 manifests; masked probes removed |
| Expression | TCGA BRCA RNA-seq, RSEM log2(x+1) (Xena HiSeqV2) | genes kept if expressed (log2 >= 1) in >= 20% of tumours |
| Subtype | PAM50 RNA-seq call (Xena BRCA clinical matrix) | Basal vs LumA+LumB for H3 |

Unit of analysis: one primary tumour (01A sample) per donor with all three assays. Tumours lacking any assay are excluded and counted.

## 3. Feasibility gate (checked before any outcome analysis)

Gate: >= 50 tumours with all three assays and >= 10 PAM50 Basal; if fewer than 10 Basal, H3 is descriptive only; if fewer than 50 tumours, stop.

Result of the check (sample identifiers and platform only, no values compared):
- 74 primary BRCA tumours have ATAC-seq; all 74 have RNA-seq.
- The UCSC Xena legacy 450K matrix covers only 43 of them. GDC lists HM27 beta files for 30 more; one tumour (TCGA-BH-A0B1) has no primary-tumour methylation file.
- Final cohort: **73 tumours** (43 HM450, 30 HM27). PAM50 (Xena RNA-seq call) over all 74: LumA 30, LumB 15, Basal 12, Her2 10, Normal 1, missing 6. Gate met.
- Platform and subtype are partly confounded (9 of 10 Her2-enriched tumours are on HM27; Basal 8 on HM450, 4 on HM27). This makes the platform-split analysis (Section 5) mandatory, and H3 must also be shown within HM450 only.

Data-source decision made at this step (before any outcome analysis): all methylation is taken from the GDC harmonized SeSAMe pipeline rather than the Xena legacy matrix, so both arrays share one processing pipeline.

## 4. Definitions (fixed now)

- **Gene set:** protein-coding GENCODE v36 basic genes; TSS = most 5' TSS of the gene's basic transcripts. Genes with several TSS > 1 kb apart: sensitivity analysis using the TSS of the most expressed transcript is not possible with gene-level RSEM, so the primary TSS is used throughout.
- **Promoter window:** TSS -1000 to +500 bp.
- **Promoter methylation:** mean beta of HM450 probes in the window; gene kept if >= 2 probes.
- **Promoter accessibility:** maximum normalized ATAC score across BRCA peaks overlapping the window; genes with no overlapping peak in any tumour are kept and assigned the cohort minimum (closed) only for H1 descriptive counts, and excluded from H2 modelling.
- **CpG island promoter:** window overlaps an HM450 manifest CpG-island annotation (HCP); otherwise LCP.
- **Variable genes for modelling:** SD of promoter beta >= 0.05 AND SD of accessibility above the cohort median SD AND expression SD >= 0.5.

## 5. Hypotheses and analyses

**H1 (descriptive, primary):** promoter methylation and accessibility are negatively correlated across tumours for most variable genes. Report the distribution of per-gene Spearman rho, the fraction with rho <= -0.3 (with tumour-bootstrap 95% CI), and counts of discordant states per gene-tumour pair (beta >= 0.6 with accessibility in the top tertile; beta <= 0.2 with accessibility in the bottom tertile), by HCP/LCP.

**H2 (primary):** for each variable gene, fit three linear models of expression across tumours: M (methylation), A (accessibility), MA (both). Performance is the out-of-sample R2 from leave-one-out CV. Primary statistics: median per-gene CV-R2 for M, A, MA; median per-gene difference MA - A and MA - M. 95% CIs by resampling tumours (1000 bootstrap replicates, whole pipeline re-run, gene selection fixed). Stratified by HCP/LCP. No direction is assumed; the result is reported whichever model wins.

**H2 purity sensitivity:** repeat H2 with tumour purity as a covariate in all three models (purity from the clinical matrix or ABSOLUTE if available; otherwise a methylation-based LUMP estimate). Stromal/immune admixture can drive both methylation and accessibility, so this is required, not optional.

**H3 (secondary):** basal-silenced genes = expression log2 fold change <= -1 and FDR < 0.05 (Welch t, BH) in Basal vs LumA+LumB. For each, classify the silencing as methylation-associated (promoter beta higher in Basal, difference >= 0.2), accessibility-associated (accessibility lower in Basal, log2 difference <= -1), both, or neither. Report proportions with Wilson CIs and compare with the same classification for luminal-silenced genes. Significance of the subtype contrast by 1000 permutations of subtype labels (tumour as the unit).

**Platform split (required):** H1 and H2 are run (a) pooled with array platform as a covariate in every model and (b) separately within HM450 and HM27 tumours. Results are interpreted only if the direction of the primary H2 contrasts agrees in both platform subsets. Agreement is defined as the same sign of the median MA - A and MA - M differences.

## 6. Negative controls

1. **Technical-only model:** predict each gene's expression from technical covariates only (ATAC sequencing depth / FRiP / TSS enrichment from Corces Data S1, 450K sample mean beta, TCGA plate). Same LOOCV. MA must exceed it to be interpreted.
2. **Mismatched-promoter control:** re-run H2 with each gene assigned the methylation and accessibility of a random other gene matched on expression mean decile (100 draws). Gives the null distribution of CV-R2.
3. **Nested CV:** any feature selection or tuning (none planned for H2) would be nested inside the CV loop.

## 7. Multiplicity and reporting

Per-gene tests in H3 use BH FDR. Effect sizes with 95% CIs are the primary results; p-values are secondary. No external cohort with all three assays exists, so results are framed as internal to TCGA, with bootstrap CIs and bootstrap optimism correction of CV-R2 where models are compared. The methylation platform split (HM450 vs HM27) provides the internal replication.

## 8. Deviations

Any change after this commit is recorded in a dated, hashed amendment before the affected analysis is run.

## 9. Code and data

Code on GitHub under the MohannadBanat organization (MIT licence) at submission; all inputs are public and listed with download dates and SHA-256 checksums.
