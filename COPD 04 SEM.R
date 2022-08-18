# Effects of Physical Activity and Physical Fitness on COPD (Aug 2022)
# Part 03: Structural Equation Modelling

# Load packages ----
rm(list = ls(all = TRUE))
setwd('C:/Users/Admin/OneDrive - University of Glasgow/MPH Project/Dataset')

library(haven)
library(dplyr)
library(descr)
library(mosaic)
library(lubridate)
library(lavaan)
library(lavaanPlot)

# Load the dataset ----
ukbtemp <- readRDS('./UKB_dataset/ukb_finalcohort_analysis.rds')
glimpse(ukbtemp)

# Outcome ----

# Use survival time from INITIAL VISIT (rather than accelerometer wear)
ukbtemp$time <- ukbtemp$time0
ukbtemp$time0 <- NULL

# Event-free survival time in months
favstats(ukbtemp$time) # All participants
favstats(ukbtemp[ukbtemp$event == 1, c("time")]) # COPD only
# Median survival: 12 yrs (whole cohort), 8 yrs (COPD)

# Among early COPD patients, ER visits and hospitalization
# increased by about 50% after 5 years
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5325091/


## Create binary variable for COPD at 8 years ----
ukbtemp <- ukbtemp %>%
  mutate(copd_yr0 = interval(date0, copd_date) / dyears(1)) %>%
  mutate(copd8 = case_when(copd_yr0 <= 8 ~ 1,
                           TRUE ~ 0)) %>%
  filter(!(time < 96 & copd8 == 0)) # Exclude LTFU: Ffup < 8 yrs and no event
freq(ukbtemp$copd8, plot = FALSE) # 179 incident COPD at 8 years

# Confirm
ukbtemp %>%
  filter(copd8 == 1) %>%
  select(id, time, copd_yr0, date0, copd_date) %>%
  head()
favstats(ukbtemp[ukbtemp$copd8 == 1, c("copd_yr0")]) # COPD dx within 8 yrs
favstats(ukbtemp[ukbtemp$copd8 == 0, c("time")]) # Min ffup 8 yrs if no COPD


# Data preparation
## RESTRICTION: No tuberculosis history ----
# Model would not converge with tuberculosis
freq(ukbtemp$tuberculosis)
ukbtemp <- ukbtemp %>% filter(tuberculosis == 0) # 69 removed

## Rescale total energy intake ----
# Warning due to variance being too large!
ukbtemp <- ukbtemp %>% mutate(energy_scaled = energy_kcal / 100)
favstats(ukbtemp$energy_scaled)

## Make BMI an ordered factor variable ----
ukbtemp$bmi_cat <- ordered(ukbtemp$bmi_cat)
ukbtemp <- ukbtemp %>% mutate(
  bmi_over = case_when(bmi_cat == "Normal" ~ 0,
                       TRUE ~ 1)
)

# Log transformations NOT used due to small variances and robust estimator


# Landmark analysis ----
ukb_allevents <- ukbtemp

## Exclude COPD diagnosis within 2 years from INITIAL VISIT
ukb <- ukb_allevents %>%
  filter(is.na(copd_yr0) | copd_yr0 >= 2)

head(ukb[!is.na(ukb$copd_yr0), c("copd_yr0", "date0", "copd_date")]) # Confirm
nrow(ukb_allevents) - nrow(ukb) # 24 excluded: Developed COPD within 2 years
nrow(ukb) # 23,519 left

rm(ukbtemp)


# SEM -----
model.struc <- '
PhysicalFitness =~ grip_ave + asmm + bmi
BehaviouralFactors =~ smoking + alcohol + energy_scaled
SociodemographicFactors =~ age + sex + dep

BehaviouralFactors ~ SociodemographicFactors
copd8 ~ PhysicalFitness + fev1 + multimorb + BehaviouralFactors + SociodemographicFactors + pollution
PhysicalFitness ~ fev1 + multimorb + BehaviouralFactors + SociodemographicFactors + pollution + spouse

bmi ~~ grip_ave + asmm + energy_scaled
smoking ~~ alcohol
'

fit.exp <- sem(model = c(model.struc),
               data = ukb, # change for sensitivity analysis
               ordered = c("copd8", "smoking", "alcohol", "sex",
                           "multimorb", "tuberculosis", "spouse",
                           "dep_cat", "bmi_over"),
               link = "probit", estimator = "WLSMV")
summary(fit.exp, fit.measures = TRUE, standardized = TRUE)


nodelabels = list(copd8 = "COPD", fev1 = "FEV1",
                  energy_scaled = "total energy intake",
                  dep = "deprivation score",
                  grip_ave = "grip strength",
                  bmi = "BMI",
                  asmm = "ASMM",
                  multimorb = "multimorbidity")


lavaanPlot(model = fit.exp, coefs = TRUE, stand = TRUE, covs = FALSE,
           node_options = list(shape = "box", fontname = "Helvetica"),
           edge_options = list(color = "grey"),
           graph_options = list(layout = "circo"),
           digits = 3, stars = "regress",
           labels = nodelabels)

SEMparstdexp <- as_tibble(standardizedSolution(fit.exp,
                                               type = "std.all",
                                               ci = TRUE, level = 0.95,
                                               se = TRUE, pvalue = TRUE,
                                               output = "data.frame"))
print(SEMparstdexp)

# Modification indices ----
# modindices(fit.exp, sort = TRUE)