# Effects of Physical Activity and Physical Fitness on COPD (Aug 2022)
# Part 03: Survival Analysis

# Load packages ----
rm(list = ls(all = TRUE))
setwd('C:/Users/Admin/OneDrive - University of Glasgow/MPH Project/Dataset')

library(haven)
library(dplyr)
library(descr)
library(mosaic)
library(ggplot2)
library(lubridate)
library(survival)
library(survminer) # Plot Cox model
library(survRM2) # Restricted mean survival time
library(gtsummary) # For tidying output
library(lmtest) # For likelihood ratio test
library(mgcv) # GAM with penalized splines
# library(gratia) # Plot mgcv GAMs
library(forestmodel) # Forest plot of HRs
library(AF) # Attributable fraction

# Load the dataset ----
ukbtemp <- readRDS('./UKB_dataset/ukb_finalcohort_analysis.rds')
glimpse(ukbtemp)


# Landmark analysis -----
favstats(ukbtemp$time) # Negative follow-up time

## Still need to exclude those with COPD at baseline (accelerometer wear)
ukb_allevents <- ukbtemp %>%
  filter(is.na(copd_date) | !(copd_date <= accel_date))

favstats(ukb_allevents$time) # All follow-up time positive
nrow(ukbtemp) - nrow(ukb_allevents) # 98 with COPD at end of accelerometer wear

## Exclude COPD diagnosis within 2 years from end of accelerometer wear time
ukb <- ukb_allevents %>%
  mutate(copd_yr = interval(accel_date, copd_date) / dyears(1)) %>%
  filter(is.na(copd_yr) | copd_yr >= 2)

head(ukb[!is.na(ukb$copd_yr), c("copd_yr", "accel_date", "copd_date")]) # Confirm
nrow(ukb_allevents) - nrow(ukb) # 63 excluded: Developed COPD within 2 years
nrow(ukb) # 23,585 left

# Exploratory data analysis: non-linearity ----
## GAM with splines ----
## ps defaults (Perperoglou 2019):
## cubic spline basis, 10 knots, third order difference penalty

cox.gam <- gam(time ~ s(accel_ave, bs = "ps") + # PA as continuous
                 s(grip_ave, bs = "ps") +
                 s(asmm_log, bs = "ps") +
                 bmi_cat +
                 s(age, bs = "ps") +
                 s(energy_kcal, bs = "ps") +
                 s(pollution_log, bs = "ps") +
                 s(fev1, bs = "ps") +
                 sex + smoking + alcohol + dep_cat + multimorb,
               data = ukb,
               family = cox.ph(),
               weights = event, # censoring info (1/0)
               method = "REML")
summary(cox.gam)

# Plot HRs ----
plot(cox.gam, all.terms = FALSE,
     trans = exp, # plot HR
     pages = 1, scale = 0,
     col = "red",
     shade = TRUE, shade.col = "gray90",
     main = "Partial effects: HRs by predictor")

# EDF ~ 1: accel_ave, grip_ave, asmm_log, age, energy_kcal
# EDF > 1: fev1 (p<0.01), pollution_log; age and fev1 smooths significantly predict the response

# p: accel_ave, asmm

# Exponentiate to reverse log
# https://stats.stackexchange.com/a/557860

# Alternative using gratia
# draw(cox.gam,
#     fun = exp, # plot HR
#     scales = "fixed",
#     pch = 1, cex = 1) 
# Takes a few minutes to plot...

rm(cox.gam)


# Create survival object -----
survobj <- Surv(ukb$time, ukb$event)
print(survobj[1:20, ])


# Kaplan-Meier estimator -----
## Between groups: PA level
plot(survfit(survobj ~ accel_cat, data = ukb),
     xlab = "Time in months", ylab = "COPD-free survival",
     main = "Kaplan-Meier")

ggsurvplot(
  fit = survfit(survobj ~ accel_cat, data = ukb),
  xlab = "Time in months", ylab = "Probability of COPD-free survival",
  title = "Kaplan-Meier Survival Curves",
  conf.int = FALSE, pval = TRUE, pval.method = TRUE,
  pval.coord = c(25,0.955), pval.method.coord = c(25,0.96),
  legend.title = "Physical activity level",
  legend.labs = c("Low", "Moderate", "High"),
  censor = FALSE,
  xlim = c(18,100), # x-axis start
  ylim = c(0.95,1), # y-axis starts at 95% probability
  risk.table = TRUE, tables.height = 0.2, tables.theme = theme_cleantable()
)


## Log-rank test -----
sd <- survdiff(survobj ~ accel_cat, data = ukb)
print(sd)
1 - pchisq(sd$chisq, length(sd$n) - 1) # p-value
rm(sd)
# p < 0.01, reject H0 that (population) survival curves are equal


# Nelson-Aalen estimator -----
# Cumulative hazard over time

plot(survfit(Surv(time, event) ~ accel_cat, data = ukb),
     conf.int = FALSE,
     ylab = "Cumulative hazard of COPD", xlab = "Time in months",
     lty = c(1,2,3), col = c(6,3,1), fun = "cumhaz")
title("Nelson-Aalen Graph")
legend("topleft", c("Low PA", "Moderate PA", "High PA"),
       lty = c(1,2,3), col = c(6,3,1))

ggsurvplot(
  fit = survfit(survobj ~ accel_cat, data = ukb),
  fun = "cumhaz",
  xlab = "Time in months", ylab = "Cumulative hazard of COPD",
  title = "Nelson-Aalen Graph",
  conf.int = FALSE,
  legend.title = "Physical activity level",
  legend.labs = c("Low", "Moderate", "High"),
  censor = FALSE,
  xlim = c(18,100), # x-axis start
  ylim = c(0,0.05),
  cumevents = TRUE, risk.table = FALSE,
  tables.height = 0.2, tables.theme = theme_cleantable()
)


# Restricted mean survival time ----
# 'The restricted mean survival time is a robust and clinically interpretable
# summary measure of the survival time distribution.
# Unlike median survival time, it is estimable even under heavy censoring.'
# (survRM2 vignette)

## Create a binary variable for exposure status
ukb <- ukb %>% 
  mutate(exp_mod = case_when(accel_cat == "Moderate" ~ 1,
                             accel_cat == "Low" ~ 0),
         exp_high = case_when(accel_cat == "High" ~ 1,
                              accel_cat == "Low" ~ 0))
freq(ukb$exp_mod, plot = FALSE)
freq(ukb$exp_high, plot = FALSE)
ukbm <- ukb %>% filter(!is.na(exp_mod))
ukbh <- ukb %>% filter(!is.na(exp_high))


## Moderate vs Low PA ----
res_mod <- rmst2(time = ukbm$time, status = ukbm$event, arm = ukbm$exp_mod,
                 tau = 96) # 8-year RMST
print(res_mod, digits = 3) # 1: Moderate, 0: Low
# plot(res_mod, xlab = "Time in months", ylab = "Probability of COPD-free survival")

## High vs Low PA ----
res_high <- rmst2(time = ukbh$time, status = ukbh$event, arm = ukbh$exp_high,
                  tau = 96) # 8-year RMST
print(res_high, digits = 3)
# plot(res_high, xlab = "Time in months", ylab = "Probability of COPD-free survival")

rm(ukbm, ukbh, res_mod, res_high)


# Cox proportional hazards model ----
## Multivariable: DAG + penalized splines ----
# P-splines penalise the integrated second derivative

# Default df = 4 more reliable vs df = 0 (AIC method), cf. Govindarajulu 2009
# nterm (knots) = 2.5*df; usually 10-40 equidistant knots (Perperoglou 2019)

### Smoothing splines: fev1, pollution_log (from EDA) -----
phmodel.dag2 <- coxph(Surv(time, event) ~ accel_cat +
                        grip_ave +
                        asmm_log +
                        pspline(fev1, df = 4) +
                        age +
                        energy_kcal +
                        pspline(pollution_log, df = 4) +
                        bmi_cat + 
                        smoking + alcohol +
                        sex + dep_cat +
                        multimorb + tuberculosis,
                      data = ukb,
                      method = "breslow")
print(phmodel.dag2)
# Spent 30 df; okay since 215 events / 30 = 7.2

### Reduce complexity by eliminating splines -----
### Remove pollution_log spline
phmodel.dag3 <- coxph(Surv(time, event) ~ accel_cat +
                        grip_ave +
                        asmm_log +
                        pspline(fev1, df = 4) +
                        age +
                        energy_kcal +
                        pollution_log +
                        bmi_cat + 
                        smoking + alcohol +
                        sex + dep_cat +
                        multimorb + tuberculosis,
                      data = ukb,
                      method = "breslow")
print(phmodel.dag3)
lrtest(phmodel.dag3, phmodel.dag2) # cannot reject H0

### Remove fev1 spline -> reduces to ALL LINEAR

## Multivariable: DAG-implied adjustment (all linear) ----
phmodel.dag <- coxph(Surv(time, event) ~ accel_cat +
                       grip_ave + bmi_cat + asmm_log + # PhysicalFitness
                       smoking + alcohol + energy_kcal + # BehaviouralFactors
                       age + sex + dep_cat + #SociodemographicFactors
                       fev1 +
                       pollution_log +
                       multimorb + tuberculosis,
                     data = ukb,
                     method = "breslow")
summary(phmodel.dag) # linear only
phmodel.dag %>% tbl_regression(exp = TRUE) # tidy output

lrtest(phmodel.dag3, phmodel.dag) # cannot reject H0


## Check PH assumption ----

### Complementary log minus log plot ----
plot(survfit(coxph(Surv(time, event) ~ accel_cat +
                     grip_ave +
                     asmm_log +
                     fev1 +
                     age +
                     energy_kcal +
                     pollution_log +
                     bmi_cat + 
                     smoking + alcohol +
                     sex +
                     dep + # deprivation as continuous
                     multimorb + tuberculosis,
                   data = ukb,
                   method = "breslow")),
     fun = "cloglog",
     xlab = "Time in months",
     ylab = "ln[-ln(survival probability)]",
     main = "log-log curves by PA level",
     lty = c(1:3), col = c(2:4), xlim = c(24,100))
legend("bottomright", c("Low PA", "Moderate PA", "High PA"),
       lty = c(1:3), col = c(2:4))
# Curves approximately parallel -> PH appears reasonable


### Schoenfeld residuals ----
# Test for non-proportionality, H0: Beta coefficient constant over time

# KM time transformation
coxres.km <- cox.zph(phmodel.dag, transform = "km", global = TRUE)
print(coxres.km)
# Log time transformation
coxres.log <- cox.zph(phmodel.dag, transform = "log", global = TRUE)
print(coxres.log) # global test not significant

# Plot Schoenfeld residuals
ggcoxzph(coxres.km, font.main = 10,
         point.size = 0.7, ggtheme = theme_bw()) # KM

ggcoxzph(coxres.log, font.main = 10,
         point.size = 0.7, ggtheme = theme_bw()) # log

# Closer look at FEV1 and age
oldpar <- par(mfrow = c(2,2), mar = c(4,4,2,2) + 0.1)

plot(coxres.km[11], xlab = "Time, months", col = c("red", "black"),
     main = "FEV1: KM")
abline(h = 0, lty = 3)
abline(h = coef(phmodel.dag)[21], col = "blue")

plot(coxres.log[11], xlab = "Time, months", col = c("red", "black"),
     main = "FEV1: Log")
abline(h = 0, lty = 3)
abline(h = coef(phmodel.dag)[21], col = "blue")

plot(coxres.km[8], xlab = "Time, months", col = c("red", "black"),
     main = "Age: KM")
abline(h = 0, lty = 3)
abline(h = coef(phmodel.dag)[10], col = "blue")

plot(coxres.log[8], xlab = "Time, months", col = c("red", "black"),
     main = "Age: Log")
abline(h = 0, lty = 3)
abline(h = coef(phmodel.dag)[10], col = "blue")

par(oldpar)
# dev.off()

# Variation in Beta(t) is small relative to Beta (constant) ->
# interpretation may not change (cf. Therneau & Grambsch 2000, p. 142)

# Plot without Schoenfeld tests
oldpar <- par(mfrow = c(3,3), mar = c(4,4,2,2) + 0.05)
plot(coxres.log, col = c("red", "black"))
par(oldpar)
# dev.off()


# Hazard ratio plots -----
# PA only
forest_model(phmodel.dag, covariates = c("accel_cat")) # No splines
# forest_model(phmodel.dag2, covariates = c("accel_cat")) # With splines


# Goodness-of-fit -----
# Harrel's c-statistic (concordance)
concordance(phmodel.dag) # no splines
concordance(phmodel.dag2) # with splines


# Interaction: smoking and PA -----
# Note: Conditioned on confounders for PA but not smoking

## Multiplicative interaction -----
# LR Test H0: Full and nested models fit equally well -> use nested

# No splines
phmodel.inter <- coxph(Surv(time, event) ~ accel_cat +
                         grip_ave +
                         bmi_cat +
                         asmm_log +
                         fev1 +
                         age +
                         energy_kcal +
                         pollution_log +
                         sex + smoking + alcohol + dep_cat +
                         multimorb + tuberculosis +
                         accel_cat:smoking, # Interaction term
                       data = ukb,
                       method = "breslow")
print(phmodel.inter) # p > 0.05 each PA category
lrtest(phmodel.inter, phmodel.dag)
# p > 0.05 -> cannot reject H0 -> no multiplicative interaction


## Additive interaction ----
# RERI_HR from VanderWeele 2011, Li and Chambless 2007
# RERI < 0 means negative interaction or less than additivity (Knol 2011)
# 'the interaction is said to be negative or “sub-additive”' (VanderWeele & Knol 2014)

# RERI & 95% CI of RERI using delta method
# Adapted from Li and Chambless
# https://stackoverflow.com/q/49077880

# Change for moderate / high PA
BETA3 <- coef(phmodel.inter)[c("accel_catHigh", "smokingEver", "accel_catModerate:smokingEver")]
V3 <- vcov(phmodel.inter)[c("accel_catHigh", "smokingEver", "accel_catModerate:smokingEver"),
                          c("accel_catHigh", "smokingEver", "accel_catModerate:smokingEver")]

varC <- BETA3[1] + BETA3[2] + BETA3[3]
varD <- exp(varC) - 1
varI <- exp(BETA3[1]) + exp(BETA3[2]) - 1
varF <- exp(BETA3[1]) + exp(BETA3[2]) - 2
varH <- exp(varC) / varD

# Partial derivatives of RERI
a1 <- exp(varC) - exp(BETA3[1])
a2 <- exp(varC) - exp(BETA3[2])
a3 <- exp(varC)
varA <- matrix(c(a1, a2, a3), nrow = 3, ncol = 1, byrow = TRUE)

reri_hr <- exp(varC) - exp(BETA3[1]) - exp(BETA3[2]) + 1
reri_var <- t(varA) %*% V3 %*% varA

# p-value for RERI
varZ <- reri_hr / sqrt(reri_var)
reri_pval = 2* (1 - pnorm(abs(varZ)))

# RERI, 95% CI and p-value
print(reri_hr + c(0,-1,1)*1.96*c(sqrt(reri_var)))
print(reri_pval)


# Population attributable fraction: Low PA -----
# Create a binary exposure variable for Low vs Not Low
ukb <- ukb %>%
  mutate(accel_catlow = case_when(accel_cat == "Low" ~ 1,
                                  accel_cat == "Moderate" |
                                    accel_cat == "High" ~ 0))
ukb$accel_catlow <- as.integer(ukb$accel_catlow)
freq(ukb$accel_catlow, plot = FALSE) # 1: low, 0: not low

# Fit a Cox model for low vs not low PA
phmodel.low <- coxph(Surv(time, event) ~ accel_catlow +
                       grip_ave + bmi_cat + asmm_log + # PhysicalFitness
                       smoking + alcohol + energy_kcal + # BehaviouralFactors
                       age + sex + dep_cat + #SociodemographicFactors
                       fev1 +
                       pollution_log +
                       multimorb + tuberculosis,
                     data = ukb,
                     method = "breslow")
summary(phmodel.low) # HR for low: 1.33

## Approximate with Miettinen formula ----
# (cf. Mansournia 2018)
paf1 <- (nrow(ukb[ukb$accel_catlow == 1, ]) / nrow(ukb)) *
  (1 - (1 / 1.33)) # HR for low PA
print(round(paf1*100, 2)) # PAF in percent

## Using AF package ----
# If 'binary' error occurs: restart R and load AF package only
paf2 <- AFcoxph(object = phmodel.low,
                data = ukb,
                exposure = "accel_catlow",
                times = c(75.58, 82.06, 88.33, 99.73)) # quartiles + max time
summary(paf2)
plot(paf2)


# Sensitivity analysis 1: PA as continuous ----
### Smoothing spline: continuous PA ----
phmodel.dag.cont <- coxph(Surv(time, event) ~
                            pspline(accel_ave, df = 0) +
                            grip_ave + bmi_cat + asmm_log + # PhysicalFitness
                            smoking + alcohol + energy_kcal + # BehaviouralFactors
                            age + sex + dep_cat + #SociodemographicFactors
                            fev1 +
                            pollution_log +
                            multimorb + tuberculosis,
                          data = ukb,
                          method = "breslow")
summary(phmodel.dag.cont) # with PA spline

# Note: Nonlinear PA penalized to 0 df with AIC method

## Plot partial effects ----
oldpar <- par(mfrow = c(3,3), mar = c(4,4,2,2) + 0.1)
termplot(phmodel.dag.cont,
         terms = c(1:2,4,7:8,11:12), # plot numeric predictors
         ylim = c(-15,15),
         se = TRUE)
par(oldpar)

# Plot continuous PA only
termplot(phmodel.dag.cont, terms = 1, se = TRUE,
         main = "Partial effect for physical activity",
         xlab = "Average vector magnitude, mg",
         ylab = "Log HR", col.term = "red", col.se = "black")


## Linear only (no PA spline) ----
phmodel.dag.cont0 <- coxph(Surv(time, event) ~
                             accel_ave +
                             grip_ave + bmi_cat + asmm_log + # PhysicalFitness
                             smoking + alcohol + energy_kcal + # BehaviouralFactors
                             age + sex + dep_cat + #SociodemographicFactors
                             fev1 +
                             pollution_log +
                             multimorb + tuberculosis,
                           data = ukb,
                           method = "breslow")
summary(phmodel.dag.cont0) # no spline


# Sensitivity analysis 2: Saturated DAG -----
# Maximal adjustment set
phmodel.max <- coxph(Surv(time, event) ~ accel_cat +
                       age + sex + spouse + # includes marital status 
                       dep_cat + pollution_log +
                       bmi_cat + asmm_log + grip_ave + fev1 +
                       smoking + alcohol + energy_kcal +
                       multimorb + tuberculosis,
                     data = ukb,
                     method = "breslow")
summary(phmodel.max) # no splines

attr(ukb$age, "label") <- "Age"
attr(ukb$asmm_log, "label") <- "Log ASMM"
attr(ukb$pollution_log, "label") <- "Log NO2 air pollution"
attr(ukb$grip_ave, "label") <- "Hand grip strength (average)"
phmodel.max %>% tbl_regression(exp = TRUE) # tidy output


# Sensitivity analysis 3: No landmark period ----
phmodel.allevents <- coxph(Surv(time, event) ~ accel_cat +
                             grip_ave + bmi_cat + asmm_log + # PhysicalFitness
                             smoking + alcohol + energy_kcal + # BehaviouralFactors
                             age + sex + dep_cat + #SociodemographicFactors
                             fev1 +
                             pollution_log +
                             multimorb + tuberculosis,
                           data = ukb_allevents,
                           method = "breslow")
summary(phmodel.allevents) # linear only
phmodel.allevents %>% tbl_regression(exp = TRUE) # tidy output
# Moderate and high PA: smaller point estimates, now both significant
