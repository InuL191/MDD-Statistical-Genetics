#!/bin/bash

# Genotype input preparation for PRS analysis
# This script summarizes the preprocessing steps used to prepare
# post-imputation genotype data for PRSice-2 and PRS-CS.
#
# Main steps:
# 1. Merge chromosome-specific post-imputation genotype files
# 2. Reconstruct .fam metadata using original raw genotype records
# 3. Update variant IDs to rsIDs
# 4. Generate PLINK binary files for downstream PRS analysis

# ==================================================
# Part 1. Merge post-imputation genotype files
# ==================================================

# Batch 1
for chr in {1..22} X; do
    echo "data/postimpqc/batch1/chr${chr}.dose.info1maf001.dedup" >> batch1_merge_list.txt
done

plink2 \
    --pmerge-list batch1_merge_list.txt \
    --make-bed \
    --out results/batch1_merged_genotype

# Batch 2
for chr in {1..22} X; do
    echo "data/postimpqc/batch2/chr${chr}.dose.info1maf001.dedup" >> batch2_merge_list.txt
done

plink2 \
    --pmerge-list batch2_merge_list.txt \
    --make-bed \
    --out results/batch2_merged_genotype


# ==================================================
# Part 2. Restore missing .fam metadata
# ==================================================

# This step was performed in R in the original workflow.
# It restores FID / sex information in merged .fam files
# by joining with the original TWB raw genotype records.

# Example R scripts can be run separately:
# Rscript src/helpers/fix_batch1_fam.R
# Rscript src/helpers/fix_batch2_fam.R


# ==================================================
# Part 3. Update variant IDs to rsIDs
# ==================================================

# Batch 1
cat data/variant_mapping/batch1/chr*_varid_rsid.dedup.updatename > results/batch1_varid_rsid.updatename

plink \
    --bfile results/batch1_merged_genotype \
    --update-name results/batch1_varid_rsid.updatename 2 1 \
    --make-bed \
    --out results/batch1_merged_genotype_updated

# Batch 2
cat data/variant_mapping/batch2/chr*_varid_rsid.dedup.updatename > results/batch2_varid_rsid.updatename

plink \
    --bfile results/batch2_merged_genotype \
    --update-name results/batch2_varid_rsid.updatename 2 1 \
    --make-bed \
    --out results/batch2_merged_genotype_updated
