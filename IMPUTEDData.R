library(readxl)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(fitdistrplus)
library(triangle)
library(purrr)
library(BaylorEdPsych)
library(mvnmle)
library(mice)
library(Amelia)
library(tidyr)
library(rstanarm)
library(plotly)
library(ggplot2)
library(gridExtra)
library(MASS)

setwd("C:/Users/USER/OneDrive/Desktop/STAT504-AM/PROJECTS")
projdata <- read_excel("C:\\Users\\USER\\OneDrive\\Desktop\\STAT504-AM\\PROJECTS\\PROJECT 1\\DataProject.xlsx")
dbltm <- read.csv("C:\\Users\\USER\\OneDrive\\Desktop\\STAT504-AM\\PROJECTS\\PROJECT 1\\DoublingTime.csv")
to_factor <- c("Gender", "AgeGroup", "Other", "Smoker", "AbstinenceStatus", "HBV", "HCV", 
               "Screening", "Diabetes", "Alcohol_consumption", "Criteria", 
               "DifuseCancer", "MetastaticCancer", "DeathStatus", "CancerStages",
               "CurativeTreatment","Thrombosis", "Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5", "Treatment6" )
projdata[to_factor] <- lapply(projdata[to_factor], factor)

projdata$Smoker <- factor(projdata$Smoker, levels = c("Smoker related cancer (Groupe A)", 
                                                      "No Smoker related Cancer (Groupe NA)"))
projdata[] <- lapply(projdata, function(x) {
  if (is.factor(x) && all(c("Yes", "No") %in% levels(x))) {
    factor(x, levels = c("Yes", "No"))
  } else {
    x
  }
})

mar_vars <- c("SurvivalInDays", "HBV", "HCV", "Other", "Screening", "TumorSize", 
              "Diabetes", "Alcohol_consumption", "CancerStages", "Thrombosis", 
              "MetastaticCancer", "CurativeTreatment", "DeathStatus", 
              "Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5", "Treatment6")

# Subset the dataset to only include the MAR variables
mar_data <- projdata[, mar_vars]

# Subset only quantitative variables
quant_vars <- c("SurvivalInDays", "TumorSize")
quant_data <- mar_data[, quant_vars]
cat_vars <- setdiff(mar_vars, quant_vars)

# Impute quantitative variables using PMM
quant_imputed <- mice(quant_data, m = 5, method = "pmm", seed = 500)


completed_quant_data <- complete(quant_imputed, 1)

# Update the original dataset with imputed quantitative variables
mar_data[, quant_vars] <- completed_quant_data
projdata[, quant_vars] <- mar_data[, quant_vars]
# Subset only the CancerStages variable (multi-level categorical)
cancer_stage_data <- mar_data[, "CancerStages", drop = FALSE]


cancer_stage_data$Dummy <- 1

# Impute CancerStages using polyreg
cancer_stage_imputed <- mice(cancer_stage_data, m = 5, method = "polyreg", seed = 500)


completed_cancer_stage_data <- complete(cancer_stage_imputed, 1)

# Remove the dummy column after imputation
completed_cancer_stage_data$Dummy <- NULL

# Update the original dataset with imputed CancerStages
mar_data[, "CancerStages"] <- completed_cancer_stage_data[, "CancerStages"]
projdata[, "CancerStages"] <- mar_data[, "CancerStages"]


binary_cat_vars <- setdiff(cat_vars, "CancerStages")  
binary_cat_data <- mar_data[, binary_cat_vars]

# Impute binary categorical variables using logreg
binary_cat_imputed <- mice(binary_cat_data, m = 5, method = "logreg", seed = 500)


completed_binary_cat_data <- complete(binary_cat_imputed, 1)

# Update the original dataset with imputed binary categorical variables
mar_data[, binary_cat_vars] <- completed_binary_cat_data
projdata[, binary_cat_vars] <- mar_data[, binary_cat_vars]

# test for MCAR
mcar_data_full <- projdata[, c("AbstinenceStatus", "Criteria", "AgeGroup", "Age", "DifuseCancer"), drop = FALSE]


mcar_result_full <- LittleMCAR(mcar_data_full)

mcar_result_full
#missingness is NOT MCAR
#IMPUTATION FOR MNAR MISSING
# Step 1: 
formula <- AbstinenceStatus ~ SurvivalInDays + HBV + HCV + Other + Screening + TumorSize + 
  Diabetes + Alcohol_consumption + CancerStages + Thrombosis + MetastaticCancer + 
  CurativeTreatment + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
  Treatment4 + Treatment5 + Treatment6

model <- stan_glm(formula, 
                  data = projdata, 
                  family = binomial,
                  prior = normal(0, 1),  # Prior for coefficients
                  prior_intercept = normal(0, 5),  # Prior for intercept
                  seed = 500)

# Step 2: 
predicted_missing_probabilities <- predict(model, type = "response")
predicted_missing_probabilities[is.na(predicted_missing_probabilities)] <- 0.5

# Step 3: 
missing_indices <- is.na(projdata$AbstinenceStatus)
imputed_values <- rbinom(sum(missing_indices), 1, predicted_missing_probabilities[missing_indices])
projdata$AbstinenceStatus[missing_indices] <- factor(imputed_values, levels = c(0, 1), labels = levels(projdata$AbstinenceStatus))

# Step 4:
formula_age_group <- AgeGroup ~ SurvivalInDays + HBV + HCV + Other + Screening + TumorSize + 
  Diabetes + Alcohol_consumption + CancerStages + Thrombosis + MetastaticCancer + 
  CurativeTreatment + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
  Treatment4 + Treatment5 + Treatment6
model_age_group <- stan_glm(formula_age_group, 
                            data = projdata, 
                            family = binomial,  
                            prior = normal(0, 1), 
                            prior_intercept = normal(0, 5), 
                            seed = 500)


predicted_missing_probabilities_age_group <- predict(model_age_group, type = "response")
predicted_missing_probabilities_age_group[is.na(predicted_missing_probabilities_age_group)] <- 0.5


missing_indices_age_group <- is.na(projdata$AgeGroup)
imputed_values_age_group <- rbinom(sum(missing_indices_age_group), 1, 
                                   predicted_missing_probabilities_age_group[missing_indices_age_group])
projdata$AgeGroup[missing_indices_age_group] <- factor(imputed_values_age_group, 
                                                       levels = c(0, 1), 
                                                       labels = levels(projdata$AgeGroup))


formula_difuse_cancer <- DifuseCancer ~ SurvivalInDays + HBV + HCV + Other + Screening + TumorSize + 
  Diabetes + Alcohol_consumption + CancerStages + Thrombosis + MetastaticCancer + 
  CurativeTreatment + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
  Treatment4 + Treatment5 + Treatment6
model_difuse_cancer <- stan_glm(formula_difuse_cancer, 
                                data = projdata, 
                                family = binomial,
                                prior = normal(0, 1), 
                                prior_intercept = normal(0, 5), 
                                seed = 500)


predicted_missing_probabilities_difuse_cancer <- predict(model_difuse_cancer, type = "response")
predicted_missing_probabilities_difuse_cancer[is.na(predicted_missing_probabilities_difuse_cancer)] <- 0.5


missing_indices_difuse_cancer <- is.na(projdata$DifuseCancer)
imputed_values_difuse_cancer <- rbinom(sum(missing_indices_difuse_cancer), 1, 
                                       predicted_missing_probabilities_difuse_cancer[missing_indices_difuse_cancer])
projdata$DifuseCancer[missing_indices_difuse_cancer] <- factor(imputed_values_difuse_cancer, 
                                                               levels = c(0, 1), 
                                                               labels = levels(projdata$DifuseCancer))


formula_criteria <- Criteria ~ AbstinenceStatus + AgeGroup + Age + DifuseCancer + 
  HBV + HCV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
  CancerStages + Thrombosis + MetastaticCancer + CurativeTreatment + DeathStatus + 
  Treatment1 + Treatment2 + Treatment3 + Treatment4 + Treatment5 + Treatment6

# Fit the Bayesian model for Criteria
model_criteria <- stan_glm(formula_criteria, 
                           data = projdata, 
                           family = "binomial",  # or "multinomial" if Criteria has multiple levels
                           prior = normal(0, 1), 
                           prior_intercept = normal(0, 5), 
                           seed = 500)


predicted_missing_probabilities_criteria <- predict(model_criteria, type = "response")


predicted_missing_probabilities_criteria[is.na(predicted_missing_probabilities_criteria)] <- 0.5


missing_indices_criteria <- is.na(projdata$Criteria)
imputed_values_criteria <- rbinom(sum(missing_indices_criteria), 1, 
                                  predicted_missing_probabilities_criteria[missing_indices_criteria])


projdata$Criteria[missing_indices_criteria] <- factor(imputed_values_criteria, 
                                                      levels = levels(projdata$Criteria))






imputed <- mice(projdata[, c("AgeGroup", "Criteria", "AbstinenceStatus")], method = "logreg", m = 5, seed = 500)


completed <- complete(imputed, 1)  # Get the first imputed dataset

projdata$AgeGroup <- completed$AgeGroup
projdata$Criteria <- completed$Criteria
projdata$AbstinenceStatus <- completed$AbstinenceStatus


summary(projdata$AgeGroup)
summary(projdata$Criteria)
summary(projdata$AbstinenceStatus)

library(randomForest)

missing_indices_age <- is.na(projdata$Age)

# Fit the random forest model
rf_model <- randomForest(Age ~ AgeGroup, data = projdata, na.action = na.exclude)


predicted_age_rf <- predict(rf_model, newdata = projdata)


projdata$Age[missing_indices_age] <- predicted_age_rf[missing_indices_age]

summary(projdata$Age)

