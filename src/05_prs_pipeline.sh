#!/bin/bash

# PRS pipeline for major depressive disorder
# This script summarizes the core workflow used to construct PRS
# using PRSice-2 (C+T) and PRS-CS for Taiwan Biobank genotype data.
#
# Notes:
# 1. Placeholder paths should be replaced with local file locations.
# 2. This repository provides a simplified and documented version
#    of the analysis workflow used in the thesis.
# 3. Only the final all-ancestry GWAS-based pipeline is retained here.

# ==================================================
# Part 1. PRSice-2: Clumping and Thresholding (C+T)
# ==================================================

# Batch 1
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base data/mdd_all_ancestry_gwas.txt.gz \
    --target data/batch1_genotype \
    --binary-target T \
    --pheno results/phenotype_batch1.txt \
    --cov results/covariates_batch1.txt \
    --base-maf EAF:0.01 \
    --stat logOR \
    --beta \
    --snp SNP \
    --chr Chromosome \
    --bp Position \
    --a1 EA \
    --a2 NEA \
    --pvalue P \
    --extract data/batch1_prs_valid_snps.txt \
    --out results/prsice_batch1_all_ancestry

# Batch 2
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base data/mdd_all_ancestry_gwas.txt.gz \
    --target data/batch2_genotype \
    --binary-target T \
    --pheno results/phenotype_batch2.txt \
    --cov results/covariates_batch2.txt \
    --base-maf EAF:0.01 \
    --stat logOR \
    --beta \
    --snp SNP \
    --chr Chromosome \
    --bp Position \
    --a1 EA \
    --a2 NEA \
    --pvalue P \
    --extract data/batch2_prs_valid_snps.txt \
    --out results/prsice_batch2_all_ancestry


# ==================================================
# Part 2. PRS-CS preprocessing
# ==================================================

# Example Python preprocessing for .bim files:
# Convert chromosome labels (e.g., X/Y) into numeric values recognized by PRS-CS.
# This step can be done separately before PRS-CS execution.

# Example for Batch 1
python - <<'PY'
import pandas as pd

bim_file = "data/batch1_genotype.bim"
output_file = "data/batch1_genotype_prscs.bim"

chr_map = {"X": 23, "Y": 24}
df = pd.read_csv(bim_file, delim_whitespace=True, header=None, dtype={0: str})
df[0] = df[0].map(lambda x: chr_map.get(x, x))
df.to_csv(output_file, sep="\t", header=False, index=False)
PY

# Example for Batch 2
python - <<'PY'
import pandas as pd

bim_file = "data/batch2_genotype.bim"
output_file = "data/batch2_genotype_prscs.bim"

chr_map = {"X": 23, "Y": 24}
df = pd.read_csv(bim_file, delim_whitespace=True, header=None, dtype={0: str})
df[0] = df[0].map(lambda x: chr_map.get(x, x))
df.to_csv(output_file, sep="\t", header=False, index=False)
PY

# Reformat GWAS summary statistics for PRS-CS input
python - <<'PY'
import pandas as pd

df = pd.read_csv("data/mdd_all_ancestry_gwas.txt", delim_whitespace=True)
df.rename(columns={
    "EA": "A1",
    "NEA": "A2",
    "logOR": "BETA",
}, inplace=True)

df["A1"] = df["A1"].str.upper()
df["A2"] = df["A2"].str.upper()

df_filtered = df[["SNP", "A1", "A2", "BETA", "P"]]
df_filtered.to_csv(
    "data/mdd_all_ancestry_gwas_prscs.txt",
    sep="\t",
    index=False,
    encoding="utf-8"
)
PY


# ==================================================
# Part 3. PRS-CS calculation
# ==================================================

module load genetics/PRScs_1.1.0

# Batch 1
python /opt/appl/PRScs_v1.1.0/PRScs.py \
    --ref_dir=data/ldblk_1kg_eas \
    --bim_prefix=data/batch1_genotype_prscs \
    --sst_file=data/mdd_all_ancestry_gwas_prscs.txt \
    --n_gwas=991073 \
    --out_dir=results/prscs_batch1_all_ancestry

# Batch 2
python /opt/appl/PRScs_v1.1.0/PRScs.py \
    --ref_dir=data/ldblk_1kg_eas \
    --bim_prefix=data/batch2_genotype_prscs \
    --sst_file=data/mdd_all_ancestry_gwas_prscs.txt \
    --n_gwas=991073 \
    --out_dir=results/prscs_batch2_all_ancestry


# ==================================================
# Part 4. Merge PRS-CS posterior effect sizes
# ==================================================

cat results/prscs_batch1_all_ancestry/*chr*.txt > results/prscs_batch1_all_ancestry_merged.txt
cat results/prscs_batch2_all_ancestry/*chr*.txt > results/prscs_batch2_all_ancestry_merged.txt


# ==================================================
# Part 5. Remove duplicated SNPs
# ==================================================

awk '{if (seen[$2]++) print $2}' data/batch1_genotype.bim > results/duplicate_snp_batch1.txt
awk '{if (seen[$2]++) print $2}' data/batch2_genotype.bim > results/duplicate_snp_batch2.txt


# ==================================================
# Part 6. Generate PLINK score files
# ==================================================

awk '{print $2, $4, $6}' results/prscs_batch1_all_ancestry_merged.txt > results/prscs_batch1_score.txt
awk '{print $2, $4, $6}' results/prscs_batch2_all_ancestry_merged.txt > results/prscs_batch2_score.txt


# ==================================================
# Part 7. Compute PRS using PLINK --score
# ==================================================

# Batch 1
plink \
    --bfile data/batch1_genotype \
    --exclude results/duplicate_snp_batch1.txt \
    --score results/prscs_batch1_score.txt 1 2 3 header \
    --a1-allele results/prscs_batch1_score.txt 2 1 \
    --out results/prscs_batch1_all_ancestry

# Batch 2
plink \
    --bfile data/batch2_genotype \
    --exclude results/duplicate_snp_batch2.txt \
    --score results/prscs_batch2_score.txt 1 2 3 header \
    --a1-allele results/prscs_batch2_score.txt 2 1 \
    --out results/prscs_batch2_all_ancestry
