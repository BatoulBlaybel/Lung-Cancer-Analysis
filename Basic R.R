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
setwd("C:/Users/USER/OneDrive/Desktop/STAT504-AM/PROJECTS")
projdata <- read_excel("C:\\Users\\USER\\OneDrive\\Desktop\\STAT504-AM\\PROJECTS\\PROJECT 1\\DataProject.xlsx")
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

#check the % of missingness in each var
missingness_summary <- sapply(projdata, function(x) sum(is.na(x)) / length(x) * 100)

# Get all categorical variables
categorical_vars <- projdata_imputed[, sapply(projdata_imputed, is.factor)]

# Create a function to compute Cramér's V for all pairs of categorical variables
compute_cramers_v <- function(data) {
  var_names <- colnames(data)
  result <- data.frame(var1 = character(0), var2 = character(0), cramer_v = numeric(0))
  
  for (i in 1:(ncol(data) - 1)) {
    for (j in (i + 1):ncol(data)) {
      # Compute Cramér's V for each pair
      cramer_v_val <- assocstats(table(data[, i], data[, j]))$cramer
      result <- rbind(result, data.frame(var1 = var_names[i], var2 = var_names[j], cramer_v = cramer_v_val))
    }
  }
  
  return(result)
}


cramers_v_matrix <- compute_cramers_v(categorical_vars)


sorted_cramers_v <- cramers_v_matrix[order(-cramers_v_matrix$cramer_v), ]
print(sorted_cramers_v)
#no high correlation between variables

#check the nature of missingness
projdata$missing_SurvivalInDays <- is.na(projdata$SurvivalInDays)
projdata$missing_AgeGroup <- is.na(projdata$AgeGroup)
projdata$missing_AbstinenceStatus <- is.na(projdata$AbstinenceStatus)
projdata$missing_HBV <- is.na(projdata$HBV)
projdata$missing_HCV <- is.na(projdata$HCV)
projdata$missing_Other <- is.na(projdata$Other)
projdata$missing_Screening <- is.na(projdata$Screening)
projdata$missing_TumorSize <- is.na(projdata$TumorSize)
projdata$missing_Diabetes <- is.na(projdata$Diabetes)
projdata$missing_Alcohol_consumption <- is.na(projdata$Alcohol_consumption)
projdata$missing_CancerStages <- is.na(projdata$CancerStages)
projdata$missing_Thrombosis <- is.na(projdata$Thrombosis)
projdata$missing_Criteria <- is.na(projdata$Criteria)
projdata$missing_DifuseCancer <- is.na(projdata$DifuseCancer)
projdata$missing_MetastaticCancer <- is.na(projdata$MetastaticCancer)
projdata$missing_CurativeTreatment <- is.na(projdata$CurativeTreatment)
projdata$missing_DeathStatus <- is.na(projdata$DeathStatus)
projdata$missing_Treatment1 <- is.na(projdata$Treatment1)
projdata$missing_Treatment2 <- is.na(projdata$Treatment2)
projdata$missing_Treatment3 <- is.na(projdata$Treatment3)
projdata$missing_Treatment4 <- is.na(projdata$Treatment4)
projdata$missing_Treatment5 <- is.na(projdata$Treatment5)
projdata$missing_Treatment6 <- is.na(projdata$Treatment6)
# Multiple imputation for high missingness

imputed_data <- mice(projdata, m = 5, method = 'pmm', maxit = 50)
projdata_imputed <- complete(imputed_data)

model_survival <- glm(missing_SurvivalInDays ~ Gender + Age + AgeGroup + Smoker + AbstinenceStatus + 
                        HBV + HCV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                        Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                        CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                        Treatment5 + Treatment6 + DeathStatus, 
                      data = projdata_imputed, family = binomial)
model_agegroup <- glm(missing_AgeGroup ~ Gender +
                       SurvivalInDays + Smoker + AbstinenceStatus + 
                        HBV + HCV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                        Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                        CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                        Treatment5 + Treatment6 + DeathStatus, 
                      data = projdata_imputed, family = binomial)
model_abstinence <- glm(missing_AbstinenceStatus ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          HBV + HCV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                          Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                          CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                          Treatment5 + Treatment6 + DeathStatus, 
                        data = projdata_imputed, family = binomial)
model_hbv <- glm(missing_HBV ~ Gender +
                   SurvivalInDays + Smoker + AgeGroup + 
                   AbstinenceStatus + HCV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                   Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                   CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                   Treatment5 + Treatment6 + DeathStatus, 
                 data = projdata_imputed, family = binomial)
model_hcv <- glm(missing_HCV ~ Gender +
                   SurvivalInDays + Smoker + AgeGroup + 
                   AbstinenceStatus + HBV + Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                   Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                   CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                   Treatment5 + Treatment6 + DeathStatus, 
                 data = projdata_imputed, family = binomial)
model_other <- glm(missing_Other ~ Gender +
                     SurvivalInDays + Smoker + AgeGroup + 
                     AbstinenceStatus + HBV + HCV + Screening + TumorSize + Diabetes + Alcohol_consumption + 
                     Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                     CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                     Treatment5 + Treatment6 + DeathStatus, 
                   data = projdata_imputed, family = binomial)
model_screening <- glm(missing_Screening ~ Gender +
                         SurvivalInDays + Smoker + AgeGroup + 
                         AbstinenceStatus + HBV + HCV + Other + TumorSize + Diabetes + Alcohol_consumption + 
                         Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                         CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                         Treatment5 + Treatment6 + DeathStatus, 
                       data = projdata_imputed, family = binomial)
model_tumor <- glm(missing_TumorSize ~ Gender +
                     SurvivalInDays + Smoker + AgeGroup + 
                     AbstinenceStatus + HBV + HCV + Other + Screening + Diabetes + Alcohol_consumption + 
                     Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                     CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                     Treatment5 + Treatment6 + DeathStatus, 
                   data = projdata_imputed, family = binomial)
model_diabetes <- glm(missing_Diabetes ~ Gender +
                        SurvivalInDays + Smoker + AgeGroup + 
                        AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Alcohol_consumption + 
                        Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                        CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                        Treatment5 + Treatment6 + DeathStatus, 
                      data = projdata_imputed, family = binomial)
model_alcohol <- glm(missing_Alcohol_consumption ~ Gender +
                       SurvivalInDays + Smoker + AgeGroup + 
                       AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                       Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
                       CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                       Treatment5 + Treatment6 + DeathStatus, 
                     data = projdata_imputed, family = binomial)
model_cancer <- glm(missing_CancerStages ~ Gender +
                      SurvivalInDays + Smoker + AgeGroup + 
                      AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                      Thrombosis + Criteria + Alcohol_consumption + DifuseCancer + MetastaticCancer + 
                      CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                      Treatment5 + Treatment6 + DeathStatus, 
                    data = projdata_imputed, family = binomial)
model_thrombosis <- glm(missing_Thrombosis ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Criteria + Alcohol_consumption + DifuseCancer + MetastaticCancer + 
                          CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                          Treatment5 + Treatment6 + DeathStatus, 
                        data = projdata_imputed, family = binomial)
model_criteria <- glm(missing_Criteria ~ Gender +
                        SurvivalInDays + Smoker + AgeGroup + 
                        AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                        CancerStages + Thrombosis + Alcohol_consumption + DifuseCancer + MetastaticCancer + 
                        CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                        Treatment5 + Treatment6 + DeathStatus, 
                      data = projdata_imputed, family = binomial)
model_diffuse <- glm(missing_DifuseCancer ~ Gender +
                       SurvivalInDays + Smoker + AgeGroup + 
                       AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                       CancerStages + Thrombosis + Alcohol_consumption + Criteria + MetastaticCancer + 
                       CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                       Treatment5 + Treatment6 + DeathStatus, 
                     data = projdata_imputed, family = binomial)
model_metastatic <- glm(missing_MetastaticCancer ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                          Treatment5 + Treatment6 + DeathStatus, 
                        data = projdata_imputed, family = binomial)
model_curative <- glm(missing_CurativeTreatment ~ Gender +
                        SurvivalInDays + Smoker + AgeGroup + 
                        AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                        CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                        MetastaticCancer + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                        Treatment5 + Treatment6 + DeathStatus, 
                      data = projdata_imputed, family = binomial)
model_death <- glm(missing_DeathStatus ~ Gender +
                     SurvivalInDays + Smoker + AgeGroup + 
                     AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                     CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                     MetastaticCancer + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
                     Treatment5 + Treatment6 + CurativeTreatment, 
                   data = projdata_imputed, family = binomial)
model_treatment1 <- glm(missing_Treatment1 ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment2 + Treatment3 + Treatment4 + 
                          Treatment5 + Treatment6 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
model_treatment2 <- glm(missing_Treatment2 ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment1 + Treatment3 + Treatment4 + 
                          Treatment5 + Treatment6 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
model_treatment3 <- glm(missing_Treatment3 ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment1 + Treatment2 + Treatment4 + 
                          Treatment5 + Treatment6 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
model_treatment4 <- glm(missing_Treatment4 ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
                          Treatment5 + Treatment6 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
model_treatment5 <- glm(missing_Treatment5 ~ Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
                          Treatment4 + Treatment6 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
model_treatment6 <- glm(missing_Treatment6 ~  Gender +
                          SurvivalInDays + Smoker + AgeGroup + 
                          AbstinenceStatus + HBV + HCV + Other + Screening + TumorSize + Diabetes + 
                          CancerStages + Thrombosis + Alcohol_consumption + Criteria + DifuseCancer + 
                          MetastaticCancer + DeathStatus + Treatment1 + Treatment2 + Treatment3 + 
                          Treatment4 + Treatment5 + CurativeTreatment, 
                        data = projdata_imputed, family = binomial)
#all factors are MAR but(age, abstinencestatus, criteria, difusecancer)