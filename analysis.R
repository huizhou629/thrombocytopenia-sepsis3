library(survival)
library(survminer)
library(tidyverse)
library(MatchIt)
library(dplyr)
library(ggplot2)
library(tableone)
library(rms)
library(writexl)
# ---------Basic Information
path <- "./"
setwd(path)
cohort <- read_csv("cohort.csv")
head(cohort)
#-------Before PSM output
dput(names(cohort))
my_vars <- c(
    "age", "gender", "diabetes", "icd_chf", "icd_renal", "icd_copd", "sofa", "saps",
    "rrt_first_day", "vent_first_day", "vaso_first_day", "creatinine", "hemoglobin", "inr", "sodium", "wbc",
    "lower_respiratory_tract_infection", "genitourinary_tract_infection", "intraabdominal_infection",
    "skin_infection", "musculoskeletal_infection", "pirmary_bacteremia_infection",
    "catheter_related_blood_stream_infection", "systemic_fungal_infection",
    "biliary_tract_infection", "cns_infection", "cardiac_pericardial_infection", "other_infection"
)
cat_vars <- c(
    "gender", "diabetes", "icd_chf", "icd_renal", "icd_copd",
    "rrt_first_day", "vent_first_day", "vaso_first_day",
    "lower_respiratory_tract_infection", "genitourinary_tract_infection", "intraabdominal_infection",
    "skin_infection", "musculoskeletal_infection", "pirmary_bacteremia_infection",
    "catheter_related_blood_stream_infection", "systemic_fungal_infection", "biliary_tract_infection",
    "cardiac_pericardial_infection", "cns_infection", "other_infection"
)
tab1 <- CreateTableOne(vars = my_vars, data = cohort, factorVars = cat_vars,
                       strata = "treat", addOverall = T)
taboutput <- print(tab1, nonnormal = NULL, missing = T, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
write.csv(taboutput, "before_psm.csv")
summary(tab1)
tab1$CatTable
tab1$ContTable

## ------------Now do PSM
fml <- my_vars %>%
    paste(collapse = " + ") %>%
    sprintf("treat ~ %s", .)
psm <- matchit(as.formula(fml), data = cohort, method = "nearest", caliper = .1, replace = F, ratio = 1) 
summary(psm)
plot(psm, type = "jitter")
mchdata <- match.data(psm)
plot(summary(psm))

## ------------After PSM do tableone
check_indicator_or <- function(indicator) {
    glm_fit <- glm(as.factor(get(indicator)) ~ treat, binomial(link = "logit"), data = mchdata)
    re <- summary(glm_fit)
    re <- data.frame(re$coefficients, "OR" = round(exp(re$coefficients[, 1]), 3), exp(confint(glm_fit)))
    colnames(re)[c(4, 6, 7)] <- c("pval", "OR_95%CI_LOW", "OR_95%CI_UP")
    re$combine <- paste(round(re$OR, 2), "(", round(re$`OR_95%CI_LOW`, 2), ",", round(re$`OR_95%CI_UP`, 2), ")", sep = "")
    re$pval <- round(re$pval, 4)
    re$indicator <- indicator
    re <- re[-1, ]
    write_xlsx(re, paste(indicator, ".xlsx", sep = ""))
}

compute_indicator_duration_stats <- function(indicator, indicator_duration) {
    vars <- c(indicator_duration)
    tab3 <- CreateTableOne(vars = vars, data = mchdata %>% filter(get(indicator) == 1),
                           factorVars = NULL, strata = "treat", addOverall = T)
    taboutput <- print(tab3, nonnormal = vars, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
    write.csv(taboutput, paste(indicator_duration, ".csv", sep = ""))
}

tab2 <- CreateTableOne(vars = my_vars, data = mchdata, factorVars = cat_vars,
                       strata = "treat", addOverall = T)
taboutput <- print(tab2, nonnormal = NULL, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
write.csv(taboutput, "after_psm.csv")
summary(tab2)
tab2$CatTable
tab2$ContTable

tab2 <- CreateTableOne(
    vars = my_vars, data = mchdata %>% filter(group != "Prevalent thrombocytopenia"), factorVars = cat_vars,
    strata = "group", addOverall = T
)
taboutput <- print(tab2, nonnormal = NULL, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
write.csv(taboutput, "after_psm_group1.csv")

tab2 <- CreateTableOne(
    vars = my_vars, data = mchdata %>% filter(group != "Incident thrombocytopenia"), factorVars = cat_vars,
    strata = "group", addOverall = T
)
taboutput <- print(tab2, nonnormal = NULL, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
write.csv(taboutput, "after_psm_group2.csv")
## Use PSM to do analysis
## did_vent, vent_duration_day, used_vaso, vaso_duration_day,
## los_icu, los_hospital, is_major_bleeding, red_blood_amount, rrt
vars <- c("los_icu", "los_hospital")
tab3 <- CreateTableOne(vars = vars, data = mchdata,
                       factorVars = NULL, strata = "treat", addOverall = T)
taboutput <- print(tab3, nonnormal = vars, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(taboutput, "los_icu_hospital.csv")

#-------- vent duration day
all_indicator <- list(
    c("did_vent", "vent_duration_day"),
    c("used_vaso", "vaso_duration_day"),
    c("is_major_bleeding", "red_blood_amount")
)

for (x in all_indicator) {
    indicator <- x[1]
    indicator_duration <- x[2]
    check_indicator_or(indicator)
    compute_indicator_duration_stats(indicator, indicator_duration)
}
check_indicator_or("rrt")

tmp_vars = c("rrt", "is_major_bleeding", "did_vent", "used_vaso")
tab2 <- CreateTableOne(vars = tmp_vars, data = mchdata, factorVars = tmp_vars,
                       strata = "treat", addOverall = T)
taboutput <- print(tab2, nonnormal = NULL, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
summary(tab2)
tab2$CatTable
tab2$ContTable

tmp_vars = c("mort_28", "mort_icu", "mort_hospital")
tab2 <- CreateTableOne(vars = tmp_vars, data = cohort, factorVars = tmp_vars,
                       strata = "group", addOverall = T)
taboutput <- print(tab2, nonnormal = NULL, exact = NULL, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
summary(tab2)
tab2$CatTable
tab2$ContTable

# -- Cox Model
cohort$gender <- factor(cohort$gender, levels = c("F", "M"))
cohort$group <- factor(cohort$group, order = F,
                          levels = c("Incident thrombocytopenia", 
                                     "Prevalent thrombocytopenia", "No thrombocytopenia"))
cohort$subgroup <- factor(cohort$subgroup, order = F, 
                          levels = c("No thrombocytopenia", "Mild thrombocytopenia", 
                                     "Moderate thrombocytopenia", "Severe thrombocytopenia"))

cox_process <- function(time, status) { 
    for (checked_var in c("treat", "subgroup", "group")) {
        cox_model <- coxph(Surv(get(time), get(status)) ~ get(checked_var) + age + gender + diabetes + icd_chf + icd_copd + icd_renal + rrt_first_day + vent_first_day + vaso_first_day + sodium + wbc + hemoglobin + creatinine + inr, data = cohort)
        re <- summary(cox_model)
        re <- as.data.frame(cbind(re$coefficients,re$conf.int))
        re$combine <- paste(round(re$`exp(coef)`, 2), "(", round(re$`lower .95`,2), ",", round(re$`upper .95`, 2), ")",sep = "")
        re$`Pr(>|z|)` <- round(re$`Pr(>|z|)`, 4)
        re$indicator <- status
        re$var <- row.names(re)
        re <- re[,c(11:12,1:10)]
        write_xlsx(re, paste(status, "_", checked_var, ".xlsx", sep = ""))
    }
}

all_mort <- list(
    c("mort_28_duration", "mort_28"),
    c("mort_icu_duration", "mort_icu"),
    c("mort_hospital_duration", "mort_hospital")
)
for (x in all_mort) {
    cox_process(x[1], x[2])
}
