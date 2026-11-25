# Internal analysis plan committed on 2026-11-25 (pan-cancer extension)

**Working title:** Promoter Chromatin Accessibility Explains Cross-Tumor Gene Expression Better Than Promoter DNA Methylation Across Cancer Types: A Pre-Specified Pan-Cancer Test of a Breast Cancer Finding in TCGA

**Study type:** secondary, observational, integrative multi-omic analysis of public tumor data (no new samples, no wet lab). The title, abstract and introduction will state this.

**Target journal:** chosen before submission; it does not affect this plan.

## 0. Status at commit — what is already known

- **Breast cancer (BRCA, 73 tumors):** fully analyzed under the internal plan committed on 2026-11-16 (SHA-256 b2b7079a…8d2d) and its deviations log; code and results archived (Zenodo version 1.0.1, https://doi.org/10.5281/zenodo.24294127). Its results are known to the authors: median per-gene leave-one-out CV-R² 0.005 for promoter methylation (M), 0.196 for promoter accessibility (A), 0.219 for both (MA); median MA − A −0.004 (jackknife 95% CI −0.011 to 0.002); A > M for 84.5% of genes. **Breast is therefore the hypothesis-generating cohort, not a test.**
- **The other cancer types:** only sample identifiers, array platform and RNA-seq sample lists (file headers) have been checked (Section 3). No methylation, accessibility or expression value of any non-breast tumor has been read. The type-specific ATAC-seq count files for all 23 types were downloaded on 2026-11-16 inside one archive, but only the BRCA file was ever read.

## 1. Question and hypotheses

Does the breast finding hold in other cancer types? Specifically, across tumors of the same cancer type:

- **H-A (primary):** promoter accessibility explains more cross-tumor expression variance than promoter methylation.
- **H-B (co-primary):** adding promoter methylation to promoter accessibility adds little explained variance.
- **H-C (secondary):** promoter methylation and accessibility are only weakly negatively correlated across tumors.

The breast basal-versus-luminal analysis (H3 of the 2026-11-16 plan) is breast-specific and is not repeated.

## 2. Data

| Assay | Source | Processing |
|---|---|---|
| ATAC-seq | Corces et al. 2018, GDC ATACseq-AWG: cancer type-specific peak sets and log2 normalized counts (`<TYPE>_log2norm.txt`) | technical replicates averaged per tumor (all aliquots of one 15-character sample averaged) |
| DNA methylation | GDC harmonized SeSAMe level-3 beta values ("Methylation Beta Value"), one primary-tumor file per sample; hg38 positions, CpG-island annotation and masks from the Zhou lab GENCODE v36 manifests; masked probes removed | see Section 4 for the probe set |
| Expression | UCSC Xena TCGA hub, RSEM log2(x+1) `TCGA.<TYPE>.sampleMap/HiSeqV2` (https://tcga.xenahubs.net; the S3 address used for breast now returns HTTP 403) | genes kept if log2 ≥ 1 in ≥ 20% of tumors of that type |
| Purity | ABSOLUTE (TCGA PanCancer Atlas) | sensitivity analysis only |

Unit of analysis: one primary tumor (sample type 01) per donor with all three assays. The 15-character TCGA sample barcode links the assays, as in the breast analysis. Metastatic and recurrent samples are excluded.

## 3. Feasibility gate (sample counts only)

Run with `P00_feasibility_counts.py` on 2026-11-25. Of 404 ATAC-seq tumor samples (404 donors, 23 types), 390 are primary; 353 primary tumors have all three assays: BRCA 73 and 280 tumors from 21 other types. GBM has no HiSeqV2 RNA-seq for its ATAC samples and drops out. Among the 280 non-breast tumors, 276 have HM450 and 4 have HM27 only (COAD 2, UCEC 2).

All-three-assay primary tumors per non-breast type: KIRP 33, PRAD 26, COAD 25 (23 HM450), LUAD 21, STAD 21, ESCA 18, LIHC 17, LUSC 16, KIRC 16, THCA 14, LGG 10, BLCA 10, ACC 9, PCPG 9, HNSC 8, TGCT 8, MESO 7, SKCM 4, UCEC 4 (2 HM450), CHOL 2, CESC 2.

**Eligible types for the per-type analysis (pre-specified):** non-breast types with ≥ 15 HM450 tumors with all three assays: KIRP, PRAD, COAD, LUAD, STAD, ESCA, LIHC, LUSC, KIRC (9 types, 191 tumors after excluding the 2 COAD HM27 tumors). Gate: ≥ 6 eligible types; met.

## 4. Definitions (fixed now; identical to breast unless stated)

- Genes: protein-coding GENCODE v36 basic; TSS = most 5' TSS; promoter window TSS −1000 to +500 bp (strand-aware).
- **Methylation probe set for non-breast types:** all unmasked HM450 probes, because the restriction to HM27 ∩ HM450 probes in breast was forced by its HM27 tumors. HM27-only tumors (4) are excluded. **Sensitivity:** the primary analysis is repeated with the breast probe set (HM27 ∩ HM450) to check that the probe set does not drive any difference from breast.
- Promoter methylation: mean beta of the window's probes; gene kept if ≥ 2 probes.
- Promoter accessibility: maximum log2 normalized ATAC score across that type's peaks overlapping the window; genes with no overlapping peak in that type are excluded for that type.
- CpG-island promoter (HCP): window overlaps a manifest CpG-island probe annotation (CGIposition = Island); otherwise LCP.
- Per type, complete cases only (breast deviation 1); variable genes: SD beta ≥ 0.05, SD accessibility > that type's median SD, SD expression ≥ 0.5.
- Models (per gene, per type): M, A, MA = linear models of expression on z-scored methylation, accessibility or both (breast deviation 4). Performance: leave-one-tumor-out CV-R² = 1 − PRESS / total sum of squares about the mean, closed form, as in the breast code (`cv_r2`).
- A type with fewer than 100 variable genes is reported but left out of the meta-analysis.

## 5. Primary and secondary analyses

**Per-type statistics** (each eligible type t): D_t = median over genes of (CV-R²_A − CV-R²_M); G_t = median over genes of (CV-R²_MA − CV-R²_A); also median CV-R² of M, A, MA and the fraction of genes with A > M. 95% CIs by delete-one-tumor jackknife (estimate ± 1.96 SE), as in breast deviation 7.

**R1 (primary, H-A):** random-effects meta-analysis (DerSimonian–Laird) of D_t across the eligible types, with jackknife SEs. H-A is supported if the pooled D is > 0 and its 95% CI excludes 0.

**R2 (consistency, H-A):** D_t > 0 in at least 8 of the 9 eligible types (one-sided sign test, p = 0.020).

**R3 (co-primary, H-B):** random-effects pooled G with 95% CI. H-B is supported if the upper 95% bound of pooled G is < 0.02 (equivalence margin: methylation adds less than 2 percentage points of explained variance).

**Secondary:**
- S1 (H-C): per type, median Spearman ρ between promoter methylation and accessibility (genes with SD beta ≥ 0.05 and a promoter peak); random-effects pooled median ρ; fraction of genes with ρ ≤ −0.3.
- S2: all non-breast tumors pooled (all 21 types, HM450 only): per gene, models with cancer-type fixed effects and predictors z-scored within type; CV-R² of M, A and MA relative to the type-only model; genes with ≥ 2 probes, a promoter peak and passing the variability filters on within-type-centered values in every contributing type; types with < 5 tumors excluded.
- S3: HCP and LCP strata of R1 and R3.
- S4: purity-adjusted R1 and R3 (ABSOLUTE purity as a covariate in all models; tumors without purity dropped).
- S5: breast re-run through the pan-cancer code (shared probe set, platform covariate) to confirm the code reproduces the published breast values; and a meta-analysis including breast, reported beside R1–R3, not replacing them.
- S6: heterogeneity: I² for R1 and R3; D_t against tumor count and mean purity (descriptive).
- S7: probe-set sensitivity (Section 4).

## 6. Negative controls (per eligible type)

1. Technical-only model: ATAC deduplicated reads (log10), TSS enrichment, FRiP and duplicate fraction (Corces Data S1), sample-mean promoter methylation; same leave-one-out CV. MA must exceed it in the type for that type's result to be interpreted.
2. Mismatched-promoter null: each gene assigned the methylation and accessibility of a random other gene in the same expression-mean decile (100 draws); median CV-R² of MA.

## 7. Multiplicity and reporting

R1 is the single primary test; R2 and R3 are pre-specified co-primary criteria and all three are reported whatever they show. Effect sizes with 95% CIs are the main results. The paper reports the breast findings as the source of the hypotheses and the 21 non-breast types as the test.

## 8. Deviations

Any change after this commit is recorded, dated and hashed, in `pancan_deviations_log.md`, stating whether it was made before or after the relevant results were seen.

## 9. Code and data

Code in the existing repository (MohannadBanat/brca-promoter-methylation-accessibility, MIT licence) as a new release; all inputs public, with download dates and SHA-256 checksums.
