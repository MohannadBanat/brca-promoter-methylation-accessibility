## Version 2.0.0: pre-specified pan-cancer test

This release adds `pancan/`: a test of the breast finding in nine further TCGA cancer types (191 primary tumors with ATAC-seq, HM450 methylation and RNA-seq), under an internal analysis plan committed on 25 November 2026 at 00:03:36 UTC (SHA-256 `c440f092afb11f7e4afb9a7204ef51550ff12d0d88ad6f5c765527fbd542fa3b`), before any non-breast data were read.

**Pre-specified criteria**

| Criterion | Result | Met |
|---|---|---|
| R1: pooled D > 0 | 0.130 (95% CI 0.053 to 0.207) | Yes |
| R2: D > 0 in ≥ 8 of 9 types | 7 of 9 (PRAD and STAD not) | **No** |
| R3: pooled G upper limit < 0.02 | −0.046 (−0.060 to −0.032) | Yes |

Also included:
- secondary analyses, negative controls and the breast re-run through the pan-cancer code (identical to the breast-stage values archived in 1.0.1);
- post hoc analyses and corrections, each logged before it was run (deviations log entries 3–5), including CV-R² re-expressed relative to the leave-one-out type mean: on the plan's scale (1 − PRESS/SST) a type-mean prediction scores 1 − (n/(n − 1))², not zero;
- Figures 1–3 and the code that builds the manuscript and its supplementary tables, with an automatic check of every number quoted in the text against the result files.

**Fix.** `scripts/00_download_data.sh` now downloads RNA-seq, clinical and legacy HM450 matrices from `https://tcga.xenahubs.net/download`, because the S3 address used in 1.0.x returns HTTP 403. The new address serves byte-identical files: SHA-256 checked against `data/SHA256SUMS_downloads.txt`.

The breast analysis code, data manifests and results are unchanged since 1.0.0.
