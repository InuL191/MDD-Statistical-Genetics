# Data preparation for phenotype and covariate files
# This script prepares phenotype and covariate files for downstream PRS analysis
# using questionnaire data, batch identifiers, and ancestry principal components.

library(dplyr)

# ==================================================
# 1. Load input data
# ==================================================

# Replace these placeholder paths with your local file locations
health_data <- read.csv("data/health_questionnaire.csv", header = TRUE, fileEncoding = "ISO-8859-1")
survey_data <- read.csv("data/survey_questionnaire.csv", header = TRUE, fileEncoding = "ISO-8859-1")
batch_data  <- read.csv("data/batch_information.csv", header = TRUE, fileEncoding = "ISO-8859-1")

# ==================================================
# 2. Merge questionnaire and batch data
# ==================================================

common_vars_health_survey <- intersect(names(health_data), names(survey_data))
survey_unique <- survey_data %>%
  select(-all_of(common_vars_health_survey[common_vars_health_survey != "Release_No"]))

merged_health_survey <- merge(health_data, survey_unique, by = "Release_No", all = TRUE)

common_vars_merged_batch <- intersect(names(merged_health_survey), names(batch_data))
batch_unique <- batch_data %>%
  select(-all_of(common_vars_merged_batch[common_vars_merged_batch != "Release_No"]))

merged_data <- merge(merged_health_survey, batch_unique, by = "Release_No", all = TRUE)

# ==================================================
# 3. Define depression phenotype across follow-up records
# ==================================================

depression_summary <- merged_data %>%
  group_by(Release_No) %>%
  summarise(
    DEPRESSION_SELF = ifelse(any(DEPRESSION_SELF == 1, na.rm = TRUE), 1, 0),
    FOLLOW = first(FOLLOW),
    .groups = "drop"
  )

analysis_base <- merged_data %>%
  select(-DEPRESSION_SELF, -FOLLOW) %>%
  distinct(Release_No, .keep_all = TRUE) %>%
  left_join(depression_summary, by = "Release_No")

# Check binary outcome distribution
table(analysis_base$DEPRESSION_SELF)

# ==================================================
# 4. Split participants into TWB batch 1 and batch 2
# ==================================================

twb1_data <- analysis_base %>%
  filter(TWB1_ID != "")

twb2_data <- analysis_base %>%
  filter(TWB2_ID != "")

# Check sample sizes
nrow(twb1_data)
nrow(twb2_data)

# Rename sample ID column to IID for downstream merging
twb1_data <- twb1_data %>%
  rename(IID = TWB1_ID)

twb2_data <- twb2_data %>%
  rename(IID = TWB2_ID)

# ==================================================
# 5. Keep variables needed for phenotype/covariate files
# ==================================================

twb1_cleaned <- twb1_data %>%
  select(IID, DEPRESSION_SELF, AGE, SEX)

twb2_cleaned <- twb2_data %>%
  select(IID, DEPRESSION_SELF, AGE, SEX)

# ==================================================
# 6. Load PCA files and merge with cleaned data
# ==================================================

pca_batch1 <- read.table("data/batch1_pca.eigenvec", header = TRUE)
pca_batch2 <- read.table("data/batch2_pca.eigenvec", header = TRUE)

raw_covariates_batch1 <- merge(twb1_cleaned, pca_batch1, by = "IID")
raw_covariates_batch2 <- merge(twb2_cleaned, pca_batch2, by = "IID")

# Reorder columns
column_order <- c(
  "FID", "IID",
  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
  "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20",
  "AGE", "SEX", "DEPRESSION_SELF"
)

raw_covariates_batch1 <- raw_covariates_batch1[, column_order]
raw_covariates_batch2 <- raw_covariates_batch2[, column_order]

# ==================================================
# 7. Create phenotype files
# ==================================================

phenotype_batch1 <- raw_covariates_batch1[, c("FID", "IID", "DEPRESSION_SELF")]
phenotype_batch2 <- raw_covariates_batch2[, c("FID", "IID", "DEPRESSION_SELF")]

# Check case counts
table(phenotype_batch1$DEPRESSION_SELF)
nrow(phenotype_batch1)
nrow(phenotype_batch2)

# Note:
# Sample size decreases after merging with PCA files,
# because only participants with available PCA results are retained.

# ==================================================
# 8. Create covariate files
# ==================================================

covariates_batch1 <- raw_covariates_batch1[, !names(raw_covariates_batch1) %in% c("DEPRESSION_SELF")]
covariates_batch2 <- raw_covariates_batch2[, !names(raw_covariates_batch2) %in% c("DEPRESSION_SELF")]

# ==================================================
# 9. Export output files
# ==================================================

write.table(phenotype_batch1, "results/phenotype_batch1.txt", quote = FALSE, row.names = FALSE)
write.table(phenotype_batch2, "results/phenotype_batch2.txt", quote = FALSE, row.names = FALSE)

write.table(covariates_batch1, "results/covariates_batch1.txt", quote = FALSE, row.names = FALSE)
write.table(covariates_batch2, "results/covariates_batch2.txt", quote = FALSE, row.names = FALSE)
