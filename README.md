# MDD-Statistical-Genetics

Risk prediction of major depressive disorder integrating genetic, epigenetic, and psychosocial factors.

---

## Project Overview

This project aims to evaluate the predictive performance of major depressive disorder (MDD) by integrating:

- Polygenic Risk Scores (PRS)
- Methylation Risk Scores (MRS)
- Psychosocial variables

Using data from the Taiwan Biobank (TWB), we construct multiple logistic regression models and compare their predictive performance.

---

## Workflow

The analysis pipeline consists of the following steps:

### 1. Phenotype and Covariate Preparation
- Merge questionnaire data (health + survey + batch)
- Define MDD phenotype (`DEPRESSION_SELF`)
- Construct covariate files including age, sex, and ancestry principal components (PCs)

Script:
src/01_phenotype_and_covariate_preparation.R

---

### 2. Genotype Preprocessing

Post-imputation genotype data (chr1–22 and X) are merged using PLINK2.

After merging, the `.fam` files may contain missing or incomplete metadata. These are restored by linking the merged genotype data back to the original TWB raw genotype records using sample IDs (IID).

Specifically:
- Missing `FID` values are replaced using original `.fam` records
- Missing or undefined `Sex` values are corrected
- The updated `.fam` file is used for downstream PRS analysis

Variant IDs are then converted to rsIDs before PRS construction.

Script:
src/03_prs_pipeline.sh

---

### 4. MRS Construction

Methylation risk scores are constructed using EWAS-derived CpG sites.

- CpG sites are selected based on p-value thresholds
- CoMeBack is used to identify co-methylated regions (CMRs)
- Representative CpG sites are selected per region
- MRS is computed based on selected CpGs

Script:
src/04_mrs_construction.R

---

### 5. Model Construction and Evaluation

Logistic regression models are built sequentially:

- Model 1: PRS only
- Model 2–5: Add psychosocial variables
- Model 6: Include interaction terms

Model performance is evaluated using metrics such as AUC.

Scripts:
src/05_prs_model_analysis.R
src/06_mrs_model_analysis.R

---

## Notes

- This repository provides a structured and reproducible version of the analysis pipeline used in the study.
- Data paths have been simplified and anonymized for reproducibility.
- Raw TWB data are not included due to access restrictions.

---

## Author

Shih-Hsiang "Shaun" Lo
Institute of Health Data Science, Department of Public Health, National Taiwan University
