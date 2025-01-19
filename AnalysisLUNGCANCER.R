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

projdata <- read_excel("C:/Users/USER/OneDrive/Desktop/STAT504-AM/PROJECTS/PROJECT 1/cleaned_projdata.xlsx")
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
#############################1- Descriptive Statistics and characteristics of patients########################################
summary_qual <- projdata %>%
  select_if(is.factor) %>%  # Select only factor variables
  gather(key = "Variable", value = "Level") %>%  # Reshape the data
  count(Variable, Level, name = "Count") %>%  # Count the occurrences of each level for each variable
  group_by(Variable) %>%  # Group by the variable
  mutate(Percentage = Count / sum(Count) * 100) %>%  # Calculate the percentage
  ungroup()  # Remove grouping

summary_qual

summary_quant <- projdata %>%
  summarise(across(
    where(is.numeric) & !starts_with("ID"),  # Select numeric columns excluding 'ID'
    list(
      mean = ~mean(.),  # Handle missing values when calculating mean
      sd = ~sd(.),  # Handle missing values when calculating standard deviation
      ci_lower = ~mean(.) - 1.96 * sd(.) / sqrt(n()),  # Calculate lower CI considering missing values
      ci_upper = ~mean(.) + 1.96 * sd(.) / sqrt(n())  # Calculate upper CI considering missing values
    ),
    .names = "{.col}_{.fn}"  # Format the output names
  ))


print(summary_quant)

# Using the Freedman-Diaconis rule, Function to calculate bin width 
calculate_binwidth <- function(data) {
  IQR_value <- IQR(data)  # Removed na.rm = TRUE as no missing values
  n <- length(data)  # No need to handle missing values
  bin_width <- 2 * IQR_value * (n ^ (-1/3))  
  return(bin_width)
}

bin_width <- calculate_binwidth(projdata$TumorSize)
print(paste("Recommended binwidth:", bin_width))

# Graphs
plot_tumor_size_interactive1 <- ggplot(projdata, aes(x = TumorSize)) +
  geom_histogram(binwidth = 10.38, fill = "purple", color = "black", alpha = 0.7) +  # Removed na.rm = TRUE
  labs(title = "Histogram of Tumor Size", x = "Tumor Size", y = "Count") +
  theme_minimal()

plot_tumor_size_interactive2 <- ggplot(projdata, aes(x = "", y = TumorSize)) +
  geom_boxplot(fill = "purple", color = "black", alpha = 0.7) +  # Removed na.rm = TRUE
  labs(title = "Boxplot of Tumor Size", x = "", y = "Tumor Size") +
  theme_minimal()

interactive_plot1 <- ggplotly(plot_tumor_size_interactive1)
interactive_plot2 <- ggplotly(plot_tumor_size_interactive2)

plot_SurvivalInDays_interactive1 <- ggplot(projdata, aes(x = SurvivalInDays)) +
  geom_histogram(binwidth = 77.59, fill = "purple", color = "black", alpha = 0.7) +  # Removed na.rm = TRUE
  labs(title = "Histogram of Survival", x = "Survival", y = "Count") +
  theme_minimal()

plot1 <- ggplotly(plot_SurvivalInDays_interactive1)

plot_SurvivalInDays_interactive2 <- ggplot(projdata, aes(x = "", y = SurvivalInDays)) +
  geom_boxplot(fill = "purple", color = "black", alpha = 0.7) +  # Removed na.rm = TRUE
  labs(title = "Boxplot of Survival", x = "", y = "Survival") +
  theme_minimal()

plot2 <- ggplotly(plot_SurvivalInDays_interactive2)

# Count number of deaths
deaths_nb <- sum(projdata$DeathStatus == "Death")

# Total population size 
population_size <- length(unique(projdata$ID))

# Calculate mortality rate
mortality_rate <- deaths_nb / population_size
mortality_rate_percent <- mortality_rate * 100

# Calculate 95% Confidence Interval using normal approximation
z_score <- 1.96 


ci_lower <- mortality_rate - z_score * sqrt((mortality_rate * (1 - mortality_rate)) / population_size)
ci_upper <- mortality_rate + z_score * sqrt((mortality_rate * (1 - mortality_rate)) / population_size)


ci_lower_percent <- ci_lower * 100
ci_upper_percent <- ci_upper * 100


cat("Mortality Rate:", round(mortality_rate_percent, 2), "%\n")
cat("95% Confidence Interval: [", round(ci_lower_percent, 2), "%, ", round(ci_upper_percent, 2), "%]\n")


# Count the number of screened individuals
screened_patient <- sum(projdata$Screening == "Screened")

# Count the number of unscreened individuals
unscreened_patient <- sum(projdata$Screening == "Unscreened")


population_size <- length(unique(projdata$ID))

# Calculate the rate of screened cancer individuals
screened_rate <- screened_patient / population_size
screened_rate_percent <- screened_rate * 100

# Calculate the rate of unscreened cancer individuals
unscreened_rate <- unscreened_patient / population_size
unscreened_rate_percent <- unscreened_rate * 100


cat("Rate of Screened Cancer:", round(screened_rate_percent, 2), "%\n")
cat("Rate of Unscreened Cancer:", round(unscreened_rate_percent, 2), "%\n")


calculate_rr_table <- function(variable_name, dataset) {
 
  levels <- unique(dataset[[variable_name]])
  

  rr_results <- data.frame(Variable = character(), 
                           Level = character(),
                           Risk = numeric(),
                           RR = numeric(),
                           stringsAsFactors = FALSE)
  
 
  if (length(levels) > 2) {
    # Set "A" as the reference level (if present)
    reference_level <- "A"
    if (!reference_level %in% levels) {
      stop(paste("Reference level", reference_level, "is not in the levels of", variable_name))
    }
    
    
    group_ref <- dataset[dataset[[variable_name]] == reference_level, ]
    risk_ref <- mean(group_ref$DeathStatus == "Death")
    
    
    for (level in levels) {
      
      group <- dataset[dataset[[variable_name]] == level, ]
      
      
      risk_group <- mean(group$DeathStatus == "Death")
      
            if (level == reference_level) {
        rr <- 1  # Reference level has RR of 1
      } else {
        rr <- risk_group / risk_ref  # Calculate RR relative to reference level
      }
      

      rr_results <- rbind(rr_results, data.frame(Variable = variable_name, Level = level, Risk = risk_group, RR = rr))
    }
  } else {
    
    for (level in levels) {
     
      group <- dataset[dataset[[variable_name]] == level, ]
      
      
      risk_group <- mean(group$DeathStatus == "Death")
      
      
      other_group <- dataset[dataset[[variable_name]] != level, ]
      risk_other <- mean(other_group$DeathStatus == "Death")
      
      
      rr <- risk_group / risk_other  # Compare the risk between the two groups
      
     
      rr_results <- rbind(rr_results, data.frame(Variable = variable_name, Level = level, Risk = risk_group, RR = rr))
    }
  }
  
  return(rr_results)
}


categorical_vars <- names(projdata)[sapply(projdata, is.factor)]  # Get column names where data is of type factor
categorical_vars <- categorical_vars[!categorical_vars %in% c("DeathStatus", "ID")]  # Remove DeathStatus and ID columns


all_rr_results <- data.frame(Variable = character(), 
                             Level = character(),
                             Risk = numeric(),
                             RR = numeric(),
                             stringsAsFactors = FALSE)


for (var in categorical_vars) {
  rr_result <- calculate_rr_table(var, projdata)
  all_rr_results <- rbind(all_rr_results, rr_result)
}


print(all_rr_results)

# Logistic Regression: Impact of TumorSize on DeathStatus
logistic_model <- glm(DeathStatus ~ TumorSize, data = projdata, family = binomial)

summary(logistic_model)

########################################################################################################################################################
# Linear Regression: Impact of TumorSize on SurvivalDays
linear_model <- lm(SurvivalInDays ~ TumorSize, data = projdata)


summary(linear_model)
###################################2- Survival Analysis|Lead time Bias###############################################
dbltm <- read.csv("C:\\Users\\USER\\OneDrive\\Desktop\\STAT504-AM\\PROJECTS\\PROJECT 1\\DoublingTime.csv")

dbltm$DT <- dbltm$DT*365.25

hist(dbltm$DT, breaks=30, col="purple", main="Tumor Doubling Time", xlab="Days", freq=FALSE)
# Gamma Distribution
gamma_fit <- fitdistr(dbltm$DT, "gamma")
ks_gamma <- ks.test(dbltm$DT, "pgamma", 1.962174675, 0.015487722)

# Weibull Distribution
weibull_fit <- fitdistr(dbltm$DT, "weibull")
ks_weibull <- ks.test(dbltm$DT, "pweibull", 1.43849737, 140.35103769)

# Exponential Distribution
exp_fit <- fitdistr(dbltm$DT, "exponential")
ks_exp <- ks.test(dbltm$DT, "pexp", 0.0078931390 )

# Lognormal Distribution
lognorm_fit <- fitdistr(dbltm$DT, "lognormal")
ks_lognorm <- ks.test(dbltm$DT, "plnorm", 4.56580355, 0.77142052)

# Compare K-S Test Results
ks_results <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  D_statistic = c(ks_gamma$statistic, ks_weibull$statistic, ks_exp$statistic, ks_lognorm$statistic),
  P_value = c(ks_gamma$p.value, ks_weibull$p.value, ks_exp$p.value, ks_lognorm$p.value)
)


print(ks_results)


log_likelihoods <- c(logLik(gamma_fit), logLik(weibull_fit), logLik(exp_fit), logLik(lognorm_fit))


aic_values <- c(AIC(gamma_fit), AIC(weibull_fit), AIC(exp_fit), AIC(lognorm_fit))


logLik_AIC_results <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  Log_Likelihood = log_likelihoods,
  AIC = aic_values
)

print(logLik_AIC_results)

data1 <- read_excel("C:/Users/USER/OneDrive/Desktop/STAT504-AM/PROJECTS/PROJECT 1/DS.xlsx")
data2 <- read_excel("C:/Users/USER/OneDrive/Desktop/STAT504-AM/PROJECTS/PROJECT 1/DNS.xlsx")

# Gamma Distribution
gamma_fit_DS <- fitdistr(data1$DS, "gamma")
gamma_sim_DS <- rgamma(length(data1$DS), 1.970091011, 0.043185846 )
ks_gamma_DS <- ks.test(data1$DS, "pgamma", 1.970091011, 0.043185846 )

# Weibull Distribution
weibull_fit_DS <- fitdistr(data1$DS, "weibull")
weibull_sim_DS <- rweibull(length(data1$DS), 1.15700373, 48.73659254)
ks_weibull_DS <- ks.test(data1$DS, "pweibull", 1.15700373, 48.73659254)

# Exponential Distribution
exp_fit_DS <- fitdistr(data1$DS, "exponential")
exp_sim_DS <- rexp(length(data1$DS), 0.021920781 )
ks_exp_DS <- ks.test(data1$DS, "pexp", 0.021920781 )

# Lognormal Distribution
lognorm_fit_DS <- fitdistr(data1$DS, "lognormal")
lognorm_sim_DS <- rlnorm(length(data1$DS), 3.54555220, 0.64966652)
ks_lognorm_DS <- ks.test(data1$DS, "plnorm", 3.54555220, 0.64966652)

# Compare K-S Test Results
ks_results_DS <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  D_statistic = c(ks_gamma_DS$statistic, ks_weibull_DS$statistic, ks_exp_DS$statistic, ks_lognorm_DS$statistic),
  P_value = c(ks_gamma_DS$p.value, ks_weibull_DS$p.value, ks_exp_DS$p.value, ks_lognorm_DS$p.value)
)


print(ks_results_DS)


log_likelihoods_DS <- c(logLik(gamma_fit_DS), logLik(weibull_fit_DS), logLik(exp_fit_DS), logLik(lognorm_fit_DS))


aic_values_DS <- c(AIC(gamma_fit_DS), AIC(weibull_fit_DS), AIC(exp_fit_DS), AIC(lognorm_fit_DS))


logLik_AIC_results_DS <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  Log_Likelihood = log_likelihoods_DS,
  AIC = aic_values_DS
)

print(logLik_AIC_results_DS)

# Gamma Distribution
gamma_fit_DNS <- fitdistr(data2$DNS, "gamma")
gamma_sim_DNS <- rgamma(length(data2$DNS), 2.200249921, 0.033368808 )
ks_gamma_DNS <- ks.test(data2$DNS, "pgamma", 2.200249921, 0.033368808 )

# Weibull Distribution
weibull_fit_DNS <- fitdistr(data2$DNS, "weibull")
weibull_sim_DNS <- rweibull(length(data2$DNS), 1.38867014, 72.88813196)
ks_weibull_DNS <- ks.test(data2$DNS, "pweibull", 1.38867014, 72.88813196)

# Exponential Distribution
exp_fit_DNS <- fitdistr(data2$DNS, "exponential")
exp_sim_DNS <- rexp(length(data2$DNS), 0.0151651290 )
ks_exp_DNS <- ks.test(data2$DNS, "pexp", 0.0151651290 )

# Lognormal Distribution
lognorm_fit_DNS <- fitdistr(data2$DNS, "lognormal")
lognorm_sim_DNS <- rlnorm(length(data2$DNS), 3.94476185, 0.70570463)
ks_lognorm_DNS <- ks.test(data2$DNS, "plnorm", 3.94476185, 0.70570463)

# Compare K-S Test Results
ks_results_DNS <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  D_statistic = c(ks_gamma_DNS$statistic, ks_weibull_DNS$statistic, ks_exp_DNS$statistic, ks_lognorm_DNS$statistic),
  P_value = c(ks_gamma_DNS$p.value, ks_weibull_DNS$p.value, ks_exp_DNS$p.value, ks_lognorm_DNS$p.value)
)


print(ks_results_DNS)


log_likelihoods_DNS <- c(logLik(gamma_fit_DNS), logLik(weibull_fit_DNS), logLik(exp_fit_DNS), logLik(lognorm_fit_DNS))


aic_values_DNS <- c(AIC(gamma_fit_DNS), AIC(weibull_fit_DNS), AIC(exp_fit_DNS), AIC(lognorm_fit_DNS))


logLik_AIC_results_DNS <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential", "Lognormal"),
  Log_Likelihood = log_likelihoods_DNS,
  AIC = aic_values_DNS
)

print(logLik_AIC_results_DNS)

set.seed(123)
M <- 1000
lognorm_fit_DT <- fitdistr(dbltm$DT, "lognormal")
DT_sim <- rlnorm(M, 4.56580355, 0.77142052)
for (i in 1:M) {
  if (DT_sim [i] > 440) {DT_sim[i]=mean(DT_sim)}
  else {DT_sim[i]=DT_sim[i]}}
summary(DT_sim) 

lognorm_fit_DS <- fitdistr(data1$DS, "lognormal")
DS_sim <- rlnorm(M, 3.54555220, 0.64966652)
for (i in 1:M) {
  if (DS_sim[i]> 100||DS_sim[i]< 10) {DS_sim[i]=mean(DS_sim)}
  else {DS_sim[i]=DS_sim[i]}}
summary(DS_sim)

lognorm_fit_DNS <- fitdistr(data2$DNS, "lognormal")
DNS_sim <- rlnorm(M, 3.94476185, 0.70570463 )
for (i in 1:M) {
  if (DNS_sim[i]> 200||DNS_sim[i]< 10) {DNS_sim[i]=mean(DNS_sim)}
  else {DNS_sim[i]=DNS_sim[i]}}
summary(DNS_sim)

LT_sim<- 3*DT_sim*((log(DNS_sim/DS_sim))/log(2))
LT_sim<- replace (LT_sim, LT_sim<0,mean(LT_sim))

summary(LT_sim)

M <- 1000
min_value <- 10
mode_value <- 12
max_value <- 14
z <- rtriangle(n = M, a = min_value, b = max_value, 
               c = mode_value)

SojournTime <- 3 * DT_sim * ((log(DNS_sim / z)) / log(2))  
SojournTime <- replace(SojournTime, SojournTime < 0, mean(SojournTime))  
SojournTime <- replace(SojournTime, SojournTime > 3650, mean(SojournTime))  
mean(SojournTime)

t <- 365
n <- 894
A2 <- NULL
B2 <- NULL
LeadTimeM22 <- NULL
SurvivalInDays <- projdata$SurvivalInDays
for (i in 1:n) {
  if (SurvivalInDays[i]<=t)
  {
    A2 [i] <-1-exp(-(1/mean(SojournTime))*SurvivalInDays[i])-(1/mean(SojournTime))*SurvivalInDays[i]*exp(-(1/mean(SojournTime))*SurvivalInDays[i])
    B2 [i] <-(1/mean(SojournTime))*(1-exp(-(1/mean(SojournTime))*SurvivalInDays[i]))
    LeadTimeM22 [i] <- A2[i]/B2[i]
  }
  else { LeadTimeM22 [i]<- (1-exp(-(1/mean(SojournTime))*t))/(1/mean(SojournTime))}
}

summary(LeadTimeM22)

CorrectedSurvival12<- SurvivalInDays-LeadTimeM22
mean(CorrectedSurvival12)

###############################2- Survival Analysis |Estimation of Mean sojourn time ################################################

gamma_fit_ST <- fitdistr(SojournTime, "gamma")
ks_gamma_ST <- ks.test(SojournTime, "pgamma",  gamma_fit_ST$estimate[1], gamma_fit_ST$estimate[2])

# Weibull Distribution
weibull_fit_ST <- fitdistr(SojournTime, "weibull")
ks_weibull_ST <- ks.test(SojournTime, "pweibull", 1.29294923, 784.86070315)

# Exponential Distribution
exp_fit_ST <- fitdistr(SojournTime, "exponential")
ks_exp_ST <- ks.test(SojournTime, "pexp", exp_fit_ST$estimate[1] )


# Compare K-S Test Results
ks_results_ST <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential"),
  D_statistic = c(ks_gamma_ST$statistic, ks_weibull_ST$statistic, ks_exp_ST$statistic),
  P_value = c(ks_gamma_ST$p.value, ks_weibull_ST$p.value, ks_exp_ST$p.value)
)


print(ks_results_ST)


log_likelihoods_ST <- c(logLik(gamma_fit_ST), logLik(weibull_fit_ST), logLik(exp_fit_ST))


aic_values_ST <- c(AIC(gamma_fit_ST), AIC(weibull_fit_ST), AIC(exp_fit_ST))


logLik_AIC_results_ST <- data.frame(
  Distribution = c("Gamma", "Weibull", "Exponential"),
  Log_Likelihood = log_likelihoods_ST,
  AIC = aic_values_ST
)

print(logLik_AIC_results_ST)

transition_rate <- 1 /mean(SojournTime)




projdata$CorrectedSurvival12 <- SurvivalInDays - LeadTimeM22
projdata$DeathStatus <- ifelse(projdata$DeathStatus == "Death", 1, 0)



factors <- setdiff(names(projdata)[sapply(projdata, is.factor)], "DeathStatus")


for (factor_name in factors) {
  
  # Kaplan-Meier Estimator: Fit the Kaplan-Meier model for each factor
  km_fit <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ as.factor(projdata[[factor_name]]), data = projdata)
  
  # Plot Kaplan-Meier survival curves
  plot_title <- paste("Kaplan-Meier Survival Curves for", factor_name)
  ggsurvplot(km_fit, data = projdata, pval = TRUE, conf.int = TRUE, 
             title = plot_title)
  
  # Log-Rank Test: Compare survival curves for different levels of the categorical factor
  log_rank_result <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ as.factor(projdata[[factor_name]]), data = projdata)
  print(paste("Log-Rank Test for", factor_name))
  print(log_rank_result)

}





surv_obj <- Surv(projdata$CorrectedSurvival12, projdata$DeathStatus)


cox_model <- coxph(
  surv_obj ~ Gender + Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Other + Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
    Treatment5 + Treatment6,
  data = projdata
)


summary(cox_model)


cox_model <- coxph(
  surv_obj ~ Gender + Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
     Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment3 + Treatment4 + 
    Treatment5 + Treatment6,
  data = projdata
)


summary(cox_model)


cox_model <- coxph(
  surv_obj ~ Gender + Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment4 + 
    Treatment5 + Treatment6,
  data = projdata
)


summary(cox_model)


cox_model <- coxph(
  surv_obj ~ Gender + Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment4 + 
    Treatment5 ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment4 + 
    Treatment5 ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2 + Treatment4  ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + AgeGroup + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2   ,
  data = projdata
)


summary(cox_model)


cox_model <- coxph(
  surv_obj ~  Age + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1 + Treatment2   ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + Smoker + AbstinenceStatus + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + Smoker  + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + Smoker  + HBV + HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~  Age + Smoker  +  HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~   Smoker  +  HCV + 
    Screening + TumorSize + Diabetes + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~   Smoker  +  HCV + 
    Screening + TumorSize + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)

cox_model <- coxph(
  surv_obj ~   Smoker  +  
    Screening + TumorSize + Alcohol_consumption + 
    Thrombosis + Criteria + CancerStages + DifuseCancer + MetastaticCancer + 
    CurativeTreatment + Treatment1    ,
  data = projdata
)


summary(cox_model)




#IN CASE THE GRAPHS OF KAPLIN MEIER ARE NOT ALL DISPLAYED (in the above for loop), WE KAN DO IT FOR EACH FACTOR SEPERATLY.


# # 1. Gender
# km_fit1 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Gender, data = projdata)
# ggsurvplot(km_fit1, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Gender")
# log_rank_result1 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Gender, data = projdata)
# print("Log-Rank Test for Gender")
# print(log_rank_result1)

# # 2. AgeGroup
# km_fit2 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ AgeGroup, data = projdata)
# ggsurvplot(km_fit2, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for AgeGroup")
# log_rank_result2 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ AgeGroup, data = projdata)
# print("Log-Rank Test for AgeGroup")
# print(log_rank_result2)

# # 3. Smoker
# km_fit3 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Smoker, data = projdata)
# ggsurvplot(km_fit3, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Smoker")
# log_rank_result3 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Smoker, data = projdata)
# print("Log-Rank Test for Smoker")
# print(log_rank_result3)

# # 4. AbstinenceStatus
# km_fit4 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ AbstinenceStatus, data = projdata)
# ggsurvplot(km_fit4, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for AbstinenceStatus")
# log_rank_result4 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ AbstinenceStatus, data = projdata)
# print("Log-Rank Test for AbstinenceStatus")
# print(log_rank_result4)

# # 5. HBV
# km_fit5 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ HBV, data = projdata)
# ggsurvplot(km_fit5, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for HBV")
# log_rank_result5 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ HBV, data = projdata)
# print("Log-Rank Test for HBV")
# print(log_rank_result5)

# # 6. HCV
# km_fit6 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ HCV, data = projdata)
# ggsurvplot(km_fit6, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for HCV")
# log_rank_result6 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ HCV, data = projdata)
# print("Log-Rank Test for HCV")
# print(log_rank_result6)

# # 7. Other
# km_fit7 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Other, data = projdata)
# ggsurvplot(km_fit7, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Other")
# log_rank_result7 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Other, data = projdata)
# print("Log-Rank Test for Other")
# print(log_rank_result7)

# # 8. Screening
# km_fit8 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Screening, data = projdata)
# ggsurvplot(km_fit8, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Screening")
# log_rank_result8 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Screening, data = projdata)
# print("Log-Rank Test for Screening")
# print(log_rank_result8)

# # 9. Diabetes
# km_fit9 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Diabetes, data = projdata)
# ggsurvplot(km_fit9, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Diabetes")
# log_rank_result9 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Diabetes, data = projdata)
# print("Log-Rank Test for Diabetes")
# print(log_rank_result9)

# # 10. Alcohol_consumption
# km_fit10 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Alcohol_consumption, data = projdata)
# ggsurvplot(km_fit10, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Alcohol_consumption")
# log_rank_result10 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Alcohol_consumption, data = projdata)
# print("Log-Rank Test for Alcohol_consumption")
# print(log_rank_result10)

# # 11. Thrombosis
# km_fit11 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Thrombosis, data = projdata)
# ggsurvplot(km_fit11, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Thrombosis")
# log_rank_result11 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Thrombosis, data = projdata)
# print("Log-Rank Test for Thrombosis")
# print(log_rank_result11)

# # 12. Criteria
# km_fit12 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Criteria, data = projdata)
# ggsurvplot(km_fit12, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Criteria")
# log_rank_result12 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Criteria, data = projdata)
# print("Log-Rank Test for Criteria")
# print(log_rank_result12)

# # 13. CancerStages
# km_fit13 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ CancerStages, data = projdata)
# ggsurvplot(km_fit13, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for CancerStages")
# log_rank_result13 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ CancerStages, data = projdata)
# print("Log-Rank Test for CancerStages")
# print(log_rank_result13)

# # 14. DifuseCancer
# km_fit14 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ DifuseCancer, data = projdata)
# ggsurvplot(km_fit14, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for DifuseCancer")
# log_rank_result14 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ DifuseCancer, data = projdata)
# print("Log-Rank Test for DifuseCancer")
# print(log_rank_result14)

# # 15. MetastaticCancer
# km_fit15 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ MetastaticCancer, data = projdata)
# ggsurvplot(km_fit15, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for MetastaticCancer")
# log_rank_result15 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ MetastaticCancer, data = projdata)
# print("Log-Rank Test for MetastaticCancer")
# print(log_rank_result15)

# # 16. CurativeTreatment
# km_fit16 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ CurativeTreatment, data = projdata)
# ggsurvplot(km_fit16, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for CurativeTreatment")
# log_rank_result16 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ CurativeTreatment, data = projdata)
# print("Log-Rank Test for CurativeTreatment")
# print(log_rank_result16)

# # 17. Treatment1
# km_fit17 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Treatment1, data = projdata)
# ggsurvplot(km_fit17, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Treatment1")
# log_rank_result17 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Treatment1, data = projdata)
# print("Log-Rank Test for Treatment1")
# print(log_rank_result17)

# # 18. Treatment2
# km_fit18 <- survfit(Surv(CorrectedSurvival12, DeathStatus) ~ Treatment2, data = projdata)
# ggsurvplot(km_fit18, data = projdata, pval = TRUE, conf.int = TRUE, 
           # title = "Kaplan-Meier Survival Curve for Treatment2")
# log_rank_result18 <- survdiff(Surv(CorrectedSurvival12, DeathStatus) ~ Treatment2, data = projdata)
# print("Log-Rank Test for Treatment2")
# print(log_rank_result18)


# # Continue similarly for other treatments...
