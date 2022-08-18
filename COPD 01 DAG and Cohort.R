# Effects of Physical Activity and Physical Fitness on COPD (Aug 2022)
# Part 01: DATA PREPARATION

# Load packages ----
rm(list = ls(all = TRUE))
setwd('C:/Users/Admin/OneDrive - University of Glasgow/MPH Project/Dataset')

library(haven)
library(dplyr)
library(descr)
library(mosaic)
library(forcats)
library(mice)
library(mice)
library(lubridate)
library(dagitty)

# Load the dataset ----
ukb <- read_dta('./UKB_dataset/ukb_subset.dta')
glimpse(ukb)


# Data DAG ----
## Specify DAG ----
copddag <- dagitty('dag {
A1AT 
smoking
alcohol
total_energy
COPD [outcome]
FEV1 
grip_strength
BMI
ASMM
age
sex
deprivation
marital_status 
multimorbidity
outdoor_air_pollution
participation
physical_activity [exposure]
tuberculosis 
A1AT -> COPD
A1AT -> FEV1
A1AT -> multimorbidity
smoking -> COPD
alcohol -> COPD
total_energy -> COPD
smoking -> FEV1
alcohol -> FEV1
total_energy -> FEV1
smoking -> grip_strength
smoking -> BMI
smoking -> ASMM
alcohol -> grip_strength
alcohol -> BMI
alcohol -> ASMM
total_energy -> grip_strength
total_energy -> BMI
total_energy -> ASMM
smoking -> multimorbidity
alcohol -> multimorbidity
total_energy -> multimorbidity
smoking -> participation
alcohol -> participation
total_energy -> participation
smoking -> physical_activity
alcohol -> physical_activity 
total_energy -> physical_activity 
smoking -> tuberculosis
alcohol -> tuberculosis
total_energy -> tuberculosis
FEV1 -> COPD
FEV1 -> grip_strength
FEV1 -> BMI
FEV1 -> ASMM
FEV1 -> physical_activity
grip_strength -> COPD
BMI -> COPD
ASMM -> COPD
grip_strength -> physical_activity
BMI -> physical_activity
ASMM -> physical_activity
age -> smoking
age -> alcohol
age -> total_energy
sex -> smoking
sex -> alcohol
sex -> total_energy
deprivation -> smoking
deprivation -> alcohol
deprivation -> total_energy
age -> COPD
sex -> COPD
deprivation -> COPD
age -> FEV1
sex -> FEV1
deprivation -> FEV1
age -> grip_strength
age -> BMI
age -> ASMM
sex -> grip_strength
sex -> BMI
sex -> ASMM
deprivation -> grip_strength
deprivation -> BMI
deprivation -> ASMM
age -> marital_status
sex -> marital_status
deprivation -> marital_status
age -> multimorbidity
sex -> multimorbidity
deprivation -> multimorbidity
age -> participation
sex -> participation
deprivation -> participation
age -> physical_activity
sex -> physical_activity
deprivation -> physical_activity
age -> tuberculosis
sex -> tuberculosis
deprivation -> tuberculosis
marital_status -> smoking
marital_status -> alcohol
marital_status -> total_energy
marital_status -> grip_strength
marital_status -> BMI
marital_status -> ASMM
marital_status -> physical_activity 
multimorbidity -> COPD 
multimorbidity -> FEV1
multimorbidity -> grip_strength
multimorbidity -> BMI
multimorbidity -> ASMM
multimorbidity -> participation
multimorbidity -> physical_activity 
multimorbidity -> tuberculosis 
outdoor_air_pollution -> COPD
outdoor_air_pollution -> FEV1
outdoor_air_pollution -> grip_strength
outdoor_air_pollution -> BMI
outdoor_air_pollution -> ASMM
outdoor_air_pollution -> physical_activity
physical_activity -> COPD 
tuberculosis -> COPD 
tuberculosis -> FEV1
tuberculosis -> grip_strength
tuberculosis -> BMI
tuberculosis -> ASMM
tuberculosis -> physical_activity
}
')

## Plot DAG ----
plot(graphLayout(copddag), abbreviate.names = FALSE)

## Adjustment sets ----
adjustmentSets(copddag, exposure = "physical_activity",
               outcome = "COPD", effect = "total")


# Missing data diagnostics -----
nonmiss <- sapply(ukb,
                  function(x)
                    sum(!is.na(x)))
missc <- sapply(ukb,
                function(x)
                  sum(is.na(x))) # Count
missp <- sapply(ukb,
                function(x)
                  mean(is.na(x))) * 100 # Percentage

# Tabulate number of missing values
missval <- as.data.frame(cbind(nonmiss, missc, missp))
rm(nonmiss, missc, missp)
colnames(missval) <- c("Non-missing", "Missing", "Missing %")
print(missval)
rm(missval)


# Convert categorical variables ------
cbind(sapply(ukb, class))

# Sex
ukb$sex <- factor(ukb$sex)
class(ukb$sex)
ukb$sex <- ukb$sex %>% recode_factor("0" = "Female", "1" = "Male")
levels(ukb$sex)

# Ethnicity
ukb$ethnicity <- factor(ukb$ethnicity)
class(ukb$ethnicity)
ukb$ethnicity <- ukb$ethnicity %>%
  recode_factor(
    "-3" = NA_character_,
    "-1" = NA_character_)
ukb <- ukb %>% mutate(
  ethnicity = fct_collapse(ethnicity,
                           "White" = c("1", "1001", "1002", "1003"),
                           "Mixed" = c("2", "2001", "2002", "2003", "2004"),
                           "Asian or Asian British" =
                             c("3", "3001", "3002", "3003", "3004", "5"),
                           "Black or Black British" =
                             c("4", "4001", "4002", "4003"),
                           "Other" = c("6")
  ))
levels(ukb$ethnicity)
freq(ukb$ethnicity, plot = FALSE)

# COPD source
ukb$copd_source <- factor(ukb$copd_source)
class(ukb$copd_source)
ukb$copd_source <- ukb$copd_source %>% recode_factor(
  "0" = "Self-reported only",
  "1" = "Hospital admission",
  "2" = "Death only",
  "11" = "Hospital primary",
  "12" = "Death primary",
  "21" = "Hospital secondary",
  "22" = "Death contributory")
levels(ukb$copd_source)

# Smoking
ukb$smoking <- factor(ukb$smoking)
class(ukb$smoking)
ukb$smoking <- ukb$smoking %>% recode_factor(
  "-3" = NA_character_,
  "0" = "Never",
  "1" = "Previous",
  "2" = "Current")
levels(ukb$smoking)

# Alcohol
ukb$alcohol <- factor(ukb$alcohol)
class(ukb$alcohol)
ukb$alcohol <- ukb$alcohol %>% recode_factor(
  "-3" = NA_character_,
  "0" = "Never",
  "1" = "Previous",
  "2" = "Current")
levels(ukb$alcohol)


# Create new variables ----

## Average grip strength -----
md.pattern(ukb[, c('grip_left', 'grip_right')]) # Check missing
ukb <- ukb %>% mutate(
  grip_ave = (grip_left + grip_right) / 2)
head(ukb$grip_ave)
# sum(!is.na(ukb$grip_ave))

## Appendicular skeletal muscle mass -----
md.pattern(ukb[, c('mass_leftarm', 'mass_rightarm',
                   'mass_leftleg', 'mass_rightleg')]) # Check missing
ukb <- ukb %>% mutate(
  asmm = mass_leftarm + mass_rightarm + mass_leftleg + mass_rightleg)
head(ukb$asmm)
# sum(!is.na(ukb$asmm))

## Age -----
# Year of birth to INITIAL VISIT
ukb <- ukb %>%
  mutate(age = year(date0) - birthyr)
favstats(ukb$age) # With values outwith 40-69 years

# Year of birth to year of accelerometer end wear
# ukb <- ukb %>%
#  mutate(age = year(accel_date) - birthyr)
# favstats(ukb$age) # min & max acceptable

## Spouse -----
# Living with husband, wife or partner
ukb <- ukb %>%
  mutate(
    spouse = case_when(
      is.na(household_0) ~ NA_real_, # missing
      if_any(.cols = starts_with("household_"),
             ~ . == -3) ~ NA_real_, # prefer not to answer
      if_any(.cols = starts_with("household_"), ~ . == 1) ~ 1, # spouse/partner
      TRUE ~ 0))

descr::freq(ukb$spouse, plot = FALSE)
ukb <- ukb %>% select(!starts_with("household_")) # Drop household columns

## End of follow-up ----
# Recode UKB centre location to nation
ukb <- ukb %>%
  mutate(
    nation = case_when(
      location %in% c("11003", "11022", "11023") ~ "Wales",
      location %in% c("11004", "11005") ~ "Scotland",
      TRUE ~ "England"))

# Add study end date
ukb <- ukb %>%
  mutate(
    date_end = case_when(
      nation == "England" ~ as.Date("2021-09-30"),
      nation == "Scotland" ~ as.Date("2021-07-21"),
      nation == "Wales" ~ as.Date("2018-02-28")))

# GENERATE ANALYTIC COHORT -----
# Confirm whether participant ID is a unique identifier
count(duplicated(ukb$id)) # no duplicates

## (1) Exclude participants without accelerometer data -----
md.pattern(ukb[, c("accel_ave", "accel_sd", "accel_good", "accel_calibrate")])
# 398,740 missing accelerometer data

ukb1 <- subset(ukb, !is.na(accel_ave))
nrow(ukb1) # 103,671 left

## (2) Exclude participants with poor accelerometer calibration ----
freq(ukb1$accel_calibrate, plot = FALSE) # 11 poor calibration

ukb2 <- subset(ukb1, accel_calibrate == 1)
nrow(ukb2) # 103,660 left

## (3) Exclude participants with poor accelerometer wear -----
freq(ukb2$accel_good, plot = FALSE) # 6,985 poor wear time

ukb3 <- subset(ukb2, accel_good == 1)
nrow(ukb3) # 96,675 left

## (3b) Exclude "implausibly high activity values" 
# Average vector magnitude > 100 mg (Ramakrishnan 2021)
print(ukb3 %>% filter(accel_ave > 100) %>% select(id, accel_ave)) # 13 rows

ukb3 <- subset(ukb3, accel_ave < 100)
nrow(ukb3) # 96,662 left

## (4) Exclude participants without grip strength measurements for both hands ----
md.pattern(ukb3[, c("grip_left", "grip_right", "grip_ave")])
# 473 missing either hand
# 186 missing both hands

ukb4 <- subset(ukb3, !is.na(ukb3$grip_left) & !is.na(ukb3$grip_right))
nrow(ukb4) # 96,189 left
md.pattern(ukb4[, c("grip_left", "grip_right", "grip_ave")])

## (5) Exclude participants with missing BMI -----
count(is.na(ukb4$bmi)) # 123 missing BMI
ukb5 <- subset(ukb4, !is.na(ukb4$bmi))
nrow(ukb5) # 96,070 left

## (6) Exclude participants with missing ASMM -----
count(is.na(ukb5$asmm)) # 1,199 missing ASMM
ukb6 <- subset(ukb5, !is.na(ukb5$asmm))
nrow(ukb6) # 94,915 left

## Missing covariates
count(is.na(ukb6$fev1)) # 23,137 missing FEV1
count(is.na(ukb6$energy1)) # 55,694 missing total energy (instance 1)

## (7) Exclude participants without acceptable FEV1 ----
ukb7 <- subset(ukb6, !is.na(ukb6$fev1))
nrow(ukb7) # 71,778 left

## (8) Exclude missing other covariates -----
ukb8 <- ukb7[
  complete.cases(ukb7[, c("age", "sex", "ethnicity",
                          "spouse", "dep", "pollution",
                          "smoking", "alcohol", "energy1")]), ]
nrow(ukb8) # 24,646 complete cases

# (9) Exclude participants with COPD at baseline ----
# Look at extremes: date of COPD diagnosis
head(freq(ukb8$copd_date))
tail(freq(ukb8$copd_date))

# Look at extremes: initial visit
head(freq(ukb8$date0))
tail(freq(ukb8$date0))

md.pattern(ukb8[, c("copd_date", "copd_source")])
# All non-missing dates have a source

# See rows with 1900-01-01 = date unknown
copd_unk <-  ukb8[ukb8$copd_date == "1900-01-01" & !is.na(ukb8$copd_date), ]
View(cbind(copd_unk[, c("id", "copd_date", "copd_source")],
           copd_unk[, c(grep("illness_", names(copd_unk)))]))
freq(copd_unk$copd_source) # All unknown dates are self-reported only
rm(copd_unk)

nrow(ukb8) # Total complete cases
count(ukb8$copd_date == "1900-01-01") # Rows with unknown COPD date
count(is.na(ukb8$copd_date)) # Rows with missing COPD date (NA = NO COPD!)

# Identify participants with unknown COPD date or COPD at baseline
ukb8 <- ukb8 %>%
  mutate(copd_baseline = case_when(copd_date == "1900-01-01" ~ 9,
                                   copd_date <= date0 ~ 1,
                                   TRUE ~ 0))
freq(ukb8$copd_baseline, plot = FALSE)

# Verify
tail(ukb8[ukb8$copd_baseline == 1, c("copd_baseline", "copd_date", "date0")])
tail(ukb8[ukb8$copd_baseline == 0, c("copd_baseline", "copd_date", "date0")])

count(ukb8$copd_baseline == 1) # 195 with COPD at baseline
count(ukb8$copd_baseline == 9) # 17 with unknown COPD date

ukb9 <- ukb8[ukb8$copd_baseline == 0, ]
nrow(ukb9) # 24,434 without COPD at baseline

# (10) Eligibility criteria ----

## Trim ages ----
## Exclude ages outwith 40-69 years at initial visit
ukb10 <- ukb9 %>%
  filter(age >= 40 & age <= 69)
favstats(ukb10$age) # 24,172 left

ukbfinal <- ukb10
rm(ukb1, ukb2, ukb3, ukb4, ukb5, ukb6, ukb7, ukb8, ukb8a, ukb8b,
   ukb9, ukb10, ukbtemp)

# Create variables for survival analysis----
## Event variable ----
# 0: no event, 1: COPD diagnosis, 2: non-COPD death
# Event has to be numeric for ggsurv() to work
ukbfinal <- ukbfinal %>%
  mutate(event = case_when(copd_date <= date_end | copd_date <= date_death ~ 1,
                           date_death <= date_end ~ 2,
                           TRUE ~ 0))
ukbfinal <- as_tibble(ukbfinal)

# Confirm
ukbfinal[ukbfinal$event == 0, c("copd_date", "date_death", "date_end", "event")] 
ukbfinal[ukbfinal$event == 1, c("copd_date", "date_death", "date_end", "event")]
ukbfinal[ukbfinal$event == 2, c("copd_date", "date_death", "date_end", "event")]

freq(ukbfinal$event, plot = FALSE) # No non-COPD deaths -> no competing risks :)


## (Event-free) survival time variable ----

# Check logical sequence of events: Death after study end?
ukbfinal %>% filter(date_death > date_end) %>% count() # none

# Assign end of follow-up (event-free)
# Censoring events: COPD, death, study end
ukbfinal <- ukbfinal %>%
  mutate(ffupend = case_when(event == 1 ~ copd_date,
                             event == 2 ~ date_death,
                             TRUE  ~ date_end))


# Calculate survival time
# From INITIAL VISIT
ukbfinal <- ukbfinal %>%
  mutate(time0 = interval(date0, ffupend) / dmonths(1)) # in months

# Confirm
print(ukbfinal[ukbfinal$event == 0, c("copd_date", "date_death", "date_end",
                                      "event", "ffupend", "time0")])
print(ukbfinal[ukbfinal$event == 1, c("copd_date", "date_death", "date_end",
                                      "event", "ffupend", "time0")])
favstats(ukbfinal$time0) # in months
freq(ukbfinal$event, plot = FALSE) # 378 incident COPD over follow-up

# From END OF WEAR
ukbfinal <- ukbfinal %>%
  mutate(time = interval(accel_date, ffupend) / dmonths(1)) # in months

# Confirm
print(ukbfinal[ukbfinal$event == 0, c("copd_date", "date_death", "date_end",
                                      "event", "ffupend", "time")])
print(ukbfinal[ukbfinal$event == 1, c("copd_date", "date_death", "date_end",
                                      "event", "ffupend", "time")])
favstats(ukbfinal$time) # in months - NOTE NEGATIVE VALUES


# Save final cohort dataset ----
saveRDS(ukbfinal, file.path(getwd(), "./UKB_dataset/ukb_finalcohort.rds"))
rm(ukb)