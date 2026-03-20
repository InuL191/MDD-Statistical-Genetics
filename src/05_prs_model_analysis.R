# PRS model analysis for major depressive disorder
# This script evaluates the association and predictive utility of
# polygenic risk scores (PRS) together with psychosocial covariates
# using logistic regression and 5-fold cross-validation.

library(dplyr)
library(pscl)
library(caret)
library(PRROC)
library(pROC)

# --------------------------------------------------
# 1. Define model specifications
# --------------------------------------------------

prs_models <- list(
  Model1 = "DEPRESSION_SELF ~ PRS_CS_ALL + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch",
  
  Model2 = "DEPRESSION_SELF ~ PRS_CS_ALL + AGE1 + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch",
  
  Model3 = "DEPRESSION_SELF ~ PRS_CS_ALL + AGE1 + SEX + MARRIAGE1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch",
  
  Model4 = "DEPRESSION_SELF ~ PRS_CS_ALL + AGE1 + SEX + MARRIAGE1 + Higher_Edu + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch",
  
  Model5 = "DEPRESSION_SELF ~ PRS_CS_ALL + AGE1 + SEX + MARRIAGE1 + Higher_Edu + INCOME_SELF + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch",
  
  Model6 = "DEPRESSION_SELF ~ PRS_CS_ALL + SEX * AGE1 + AGE1 * MARRIAGE1 + SEX * MARRIAGE1 + AGE1 * INCOME_SELF + AGE1 * Higher_Edu +
                              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                              PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              batch"
)

cleaning_steps <- list(
  Model2 = c("AGE1", "SEX"),
  Model3 = c("MARRIAGE1"),
  Model4 = c("EDUCATION"),
  Model5 = c("INCOME_SELF")
)

# --------------------------------------------------
# 2. Helper function for incremental data cleaning
# --------------------------------------------------

prepare_model_data <- function(data, model_name, cleaning_steps) {
  cleaned_data <- data
  
  if (model_name %in% names(cleaning_steps)) {
    for (var in cleaning_steps[[model_name]]) {
      
      if (var == "EDUCATION") {
        cleaned_data <- cleaned_data %>%
          filter(!is.na(EDUCATION) & EDUCATION != "" & EDUCATION != "N" & EDUCATION != "R") %>%
          mutate(Higher_Edu = ifelse(EDUCATION %in% c(6, 7), 1, 0)) %>%
          mutate(Higher_Edu = relevel(factor(Higher_Edu), ref = "1")) %>%
          select(-EDUCATION)
        
      } else if (var == "SEX") {
        cleaned_data <- cleaned_data %>%
          filter(!is.na(SEX) & SEX != "" & SEX != "N" & SEX != "R") %>%
          mutate(SEX = relevel(factor(SEX), ref = "1"))
        
      } else if (var == "MARRIAGE1") {
        cleaned_data <- cleaned_data %>%
          filter(!is.na(MARRIAGE1) & MARRIAGE1 != "" & MARRIAGE1 != "N" & MARRIAGE1 != "R") %>%
          mutate(MARRIAGE1 = relevel(factor(MARRIAGE1), ref = "2"))
        
      } else if (var == "INCOME_SELF") {
        cleaned_data <- cleaned_data %>%
          filter(!is.na(INCOME_SELF) & INCOME_SELF != "" & INCOME_SELF != "N" & INCOME_SELF != "R") %>%
          mutate(INCOME_SELF = as.numeric(INCOME_SELF))
        
      } else {
        cleaned_data <- cleaned_data %>%
          filter(!is.na(.data[[var]]) & .data[[var]] != "" & .data[[var]] != "N" & .data[[var]] != "R")
      }
    }
  }
  
  return(cleaned_data)
}

# --------------------------------------------------
# 3. Logistic regression analysis
# --------------------------------------------------

run_prs_logistic_models <- function(data, models, cleaning_steps) {
  data_cleaned <- data
  
  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    formula_str <- as.formula(models[[i]])
    
    cat("\n============================\n")
    cat("Running", model_name, "\n")
    cat("============================\n")
    
    data_cleaned <- prepare_model_data(data_cleaned, model_name, cleaning_steps)
    
    observations <- nrow(data_cleaned)
    cat("\nObservations after cleaning:", observations, "\n")
    
    glm_model <- glm(
      formula = formula_str,
      data = data_cleaned,
      family = binomial(link = "logit")
    )
    
    print(summary(glm_model))
    
    cat("\nNagelkerke R2:", pR2(glm_model)[6], "\n")
    cat("\nNull log-likelihood:", logLik(glm(update(glm_model, . ~ 1)))[1], "\n")
    cat("Fitted log-likelihood:", logLik(glm_model)[1], "\n")
    cat("\nModel degrees of freedom:", length(coef(glm_model)) - 1, "\n")
    cat("Residual degrees of freedom:", glm_model$df.residual, "\n")
  }
}

# --------------------------------------------------
# 4. Cross-validation performance evaluation
# --------------------------------------------------

run_prs_cv_evaluation <- function(data, models, cleaning_steps, k = 5) {
  data_cleaned <- data
  
  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    formula_str <- as.formula(models[[i]])
    
    cat("\n============================\n")
    cat("Evaluating", model_name, "\n")
    cat("============================\n")
    
    data_cleaned <- prepare_model_data(data_cleaned, model_name, cleaning_steps)
    
    observations <- nrow(data_cleaned)
    cat("\nObservations after cleaning:", observations, "\n")
    
    glm_model <- glm(
      formula = formula_str,
      data = data_cleaned,
      family = binomial(link = "logit")
    )
    
    cat("\nNagelkerke R2:", pR2(glm_model)[6], "\n")
    
    folds <- createFolds(data_cleaned$DEPRESSION_SELF, k = k, list = TRUE, returnTrain = TRUE)
    
    pr_auc_list <- c()
    roc_auc_list <- c()
    
    for (fold in folds) {
      train_data <- data_cleaned[fold, ]
      test_data  <- data_cleaned[-fold, ]
      
      cv_model <- glm(
        formula = formula_str,
        data = train_data,
        family = binomial(link = "logit")
      )
      
      test_pred_prob <- predict(cv_model, newdata = test_data, type = "response")
      test_labels <- test_data$DEPRESSION_SELF
      
      pr <- pr.curve(
        scores.class0 = test_pred_prob[test_labels == 1],
        scores.class1 = test_pred_prob[test_labels == 0],
        curve = TRUE
      )
      pr_auc_list <- c(pr_auc_list, pr$auc.integral)
      
      roc_obj <- roc(test_labels, test_pred_prob)
      roc_auc_list <- c(roc_auc_list, auc(roc_obj))
    }
    
    cat("\n5-fold CV ROC-AUC:", mean(roc_auc_list, na.rm = TRUE), "\n")
    cat("5-fold CV PR-AUC:", mean(pr_auc_list, na.rm = TRUE), "\n")
  }
}

# --------------------------------------------------
# 5. Run analysis
# --------------------------------------------------

run_prs_logistic_models(
  data = TWB_merged_final,
  models = prs_models,
  cleaning_steps = cleaning_steps
)

run_prs_cv_evaluation(
  data = TWB_merged_final,
  models = prs_models,
  cleaning_steps = cleaning_steps,
  k = 5
)
