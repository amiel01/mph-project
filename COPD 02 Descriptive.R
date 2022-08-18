# Effects of Physical Activity and Physical Fitness on COPD (Aug 2022)
# Part 02: Descriptive Statistics

# Load packages ----
rm(list = ls(all = TRUE))
setwd('C:/Users/Admin/OneDrive - University of Glasgow/MPH Project/Dataset')

library(haven)
library(tidyr)
library(dplyr)
library(forcats)
library(mice)
library(descr)
library(mosaic)
library(Hmisc)
library(ggplot2)
library(gridExtra) # Multiple ggplots
library(tableone)

# Load the dataset ----
ukb <- readRDS('./UKB_dataset/ukb_finalcohort.rds')
glimpse(ukb) # 24,135 rows

# Variables
# NUMERIC: PA, hand grip strength, ASMM, pollution, total energy
# CATEGORICAL: smoking, alcohol, ethnicity, BMIcat, depcat 
# BINARY: sex, multimorbidity

# NUMERIC -----

# Convert total energy intake: kJ -> kcal
# Ref (10 Jul 2022): FAO, https://www.fao.org/3/Y5022E/y5022e04.htm
ukb$energy_kcal <- round(ukb$energy1 / 4.184, digits = 3)
print(ukb[, c("energy_kcal", "energy1")])

## Summarise all -----
# sapply(select(ukb, where(is.numeric)), favstats)
allnumeric <- sapply(ukb %>%
                       select("accel_ave", "grip_ave", "asmm", "bmi",
                              "age", "fev1", "pollution", "energy_kcal"),
                     favstats)
print(t(allnumeric))
# Zero values for grip strength and total energy
# Very high total energy values
# Need to trim


## Average vector magnitude -----
ggplot(data = ukb) +
  geom_histogram(aes(x = accel_ave),
                 colour = 1,
                 binwidth = 2) +
  labs(title = "Accelerometer PA",
       y = "Frequency", x = "Average vector magnitude (milligravity)") +
  theme_light() +
  coord_cartesian(xlim = c(0,120))

qqnorm(ukb$accel_ave, main = "Average vector magnitude")
qqline(ukb$accel_ave)
# Non-normal: categorise

### Divide Accelerometer PA into tertiles ----
ukb <- ukb %>% mutate(
  accel_cat = ntile(accel_ave, 3)) %>%
  arrange(accel_cat)

ukb$accel_cat <- ukb$accel_cat %>%
  recode_factor("1" = "Low",
                "2" = "Moderate",
                "3" = "High")
freq(ukb$accel_cat, plot = FALSE)


## Hand grip strength -----
# Zero values found

### Trim grip strength using z-scores -----
ukb <- ukb %>%
  mutate(grip_zscore = (grip_ave - mean(grip_ave, na.rm = TRUE)) /
           sd(grip_ave, na.rm = TRUE))

count(ukb[abs(ukb$grip_zscore) > 3, ]) # 46 outliers
ukb <- ukb %>% filter(!abs(grip_zscore) > 3)
nrow(ukb) # 24,089 left

ggplot(data = ukb) +
  geom_histogram(aes(x = grip_ave),
                 colour = 1,
                 binwidth = 5) +
  labs(title = "Hand grip strength: After trimming",
       y = "Frequency", x = "Average hand grip strength (kg)") +
  theme_light() 

qqnorm(ukb$grip_ave, main = "Hand grip strength")
qqline(ukb$grip_ave)
# Normality reasonable


## Appendicular skeletal muscle mass -----
### Trim ASMM using z-scores -----
ukb <- ukb %>%
  mutate(asmm_zscore =  (asmm - mean(asmm)) / sd(asmm))

nrow(ukb)
count(ukb[abs(ukb$asmm_zscore) > 3, ]) # 75 outliers
ukb <- ukb %>% filter(!abs(asmm_zscore) > 3)
nrow(ukb) # 24,089 left

hasmm <-
  ggplot(data = ukb) +
  geom_histogram(aes(x = asmm),
                 colour = 1,
                 binwidth = 2) +
  labs(title = "ASMM: After trimming",
       y = "Frequency", x = "Appendicular skeletal muscle mass (kg)") +
  theme_light()
print(hasmm)

qqnorm(ukb$asmm, main = "Appendicular skeletal muscle mass:\nAfter trimming")
qqline(ukb$asmm)
# Non-normal (bimodal)

### Transform ASMM ----
ukb$asmm_log <- log(ukb$asmm)
ukb$asmm_sqrt <- sqrt(ukb$asmm)
ukb$asmm_cubert <- (ukb$asmm)^(1/3)

# Histograms
hasmm1 <-
  ggplot(data = ukb) +
  geom_histogram(aes(x = asmm_log),
                 colour = 1,
                 binwidth = 0.1) +
  labs(title = "Log ASMM",
       y = "Frequency", x = "Log of appendicular skeletal muscle mass (kg)") +
  theme_light() 

hasmm2 <-
  ggplot(data = ukb) +
  geom_histogram(aes(x = asmm_sqrt),
                 colour = 1,
                 binwidth = 0.1) +
  labs(title = "Square root ASMM",
       y = "Frequency", x = "Square root of appendicular skeletal muscle mass (kg)") +
  theme_light() 

hasmm3 <-
  ggplot(data = ukb) +
  geom_histogram(aes(x = asmm_cubert),
                 colour = 1,
                 binwidth = 0.1) +
  labs(title = "Cube root ASMM",
       y = "Frequency", x = "Cube root of appendicular skeletal muscle mass (kg)") +
  theme_light() 

grid.arrange(hasmm, hasmm1, hasmm2, hasmm3)
rm(hasmm, hasmm1, hasmm2, hasmm3)

# Boxplots
# Plotting area in 2*2 array, # of lines of margin (bottom, left, top, right)
oldpar <- par(mfrow = c(2,2), mar = c(4,4,2,2) + 0.1)

qqnorm(ukb$asmm, main = "ASMM")
qqline(ukb$asmm)

qqnorm(ukb$asmm_log, main = "Log ASMM")
qqline(ukb$asmm_log)

qqnorm(ukb$asmm_sqrt, main = "Sqrt ASMM")
qqline(ukb$asmm_sqrt)

qqnorm(ukb$asmm_cubert, main = "Cube root ASMM")
qqline(ukb$asmm_cubert)

par(oldpar) # Switch back to default graph settings
# dev.off

# Stick with log transformed ASMM for analysis


## BMI -----
ggplot(data = ukb) +
  geom_histogram(aes(x = bmi),
                 colour = 1) +
  labs(title = "Body mass index",
       y = "Frequency", x = "BMI (kg/m^2)") +
  theme_light()

qqnorm(ukb$bmi, main = "Body mass index")
qqline(ukb$bmi)
# Non-normal: Categorise BMI

### BMI: WHO/NHS Categories ------
# Ref: nhs.uk/common-health-questions/lifestyle/what-is-the-body-mass-index-bmi
# Ref: who.int/news-room/fact-sheets/detail/obesity-and-overweight
ukb <- ukb %>%
  mutate(
    bmi_cat = case_when(bmi < 25 ~ "Normal",
                        bmi < 30 ~ "Overweight",
                        TRUE ~ "Obese"),
    bmi_cat2 = case_when(bmi < 18.5 ~ "Underweight",
                         bmi < 25 ~ "Normal",
                         bmi < 30 ~ "Overweight",
                         TRUE ~ "Obese"))

# Without underweight
ukb$bmi_cat <- as.factor(ukb$bmi_cat)
freq(ukb$bmi_cat, plot = FALSE)

# With underweight (small values)
ukb$bmi_cat2 <- as.factor(ukb$bmi_cat2)
freq(ukb$bmi_cat2, plot = FALSE)


## Age -----
ggplot(data = ukb) +
  geom_histogram(aes(x = age),
                 colour = 1,
                 binwidth = 5) +
  labs(title = "Age",
       y = "Frequency", x = "Age (years)") +
  theme_light()

qqnorm(ukb$age, main = "Age")
qqline(ukb$age)
# Reasonably normal with some deviations at extremes


## FEV1 -----
ggplot(data = ukb) +
  geom_histogram(aes(x = fev1),
                 colour = 1,
                 binwidth = 0.5) +
  labs(title = "FEV1",
       y = "Frequency", x = "FEV1 (L/min)") +
  theme_light()

qqnorm(ukb$fev1, main = "FEV1")
qqline(ukb$fev1)
# Reasonably normal


## Deprivation -----
ggplot(data = ukb) +
  geom_histogram(aes(x = dep),
                 colour = 1,
                 binwidth = 1) +
  labs(title = "Deprivation score",
       y = "Frequency", x = "Townsend deprivation index)") +
  theme_light()

qqnorm(ukb$dep, main = "Deprivation score")
qqline(ukb$dep)
# Non-normal: categorise deprivation score

### Divide deprivation into deciles ----
ukb <- ukb %>% mutate(
  dep_cat = ntile(dep, 10))
freq(ukb$dep_cat, plot = FALSE)


## Air pollution -----
ggplot(data = ukb) +
  geom_histogram(aes(x = pollution),
                 colour = 1,
                 binwidth = 1) +
  labs(title = "Air pollution",
       y = "Frequency", x = "Nitrogen dioxide air pollution (micro-g/m3)") +
  theme_light()

qqnorm(ukb$pollution, main = "Air pollution")
qqline(ukb$pollution)
# Non-normal

### Log transform air pollution -----
ukb$pollution_log <- log(ukb$pollution)

ggplot(data = ukb) +
  geom_histogram(aes(x = pollution_log),
                 colour = 1,
                 binwidth = 0.1) +
  labs(title = "Log of air pollution",
       y = "Frequency", x = "Log of nitrogen dioxide air pollution (micro-g/m3)") +
  theme_light()

qqnorm(ukb$pollution_log, main = "Log of air pollution")
qqline(ukb$pollution_log)


## Total energy intake -----
# Zero values and very high values

### Trim energy using z-scores -----
ukb <- ukb %>%
  mutate(energy_zscore = (energy_kcal - mean(energy_kcal)) / sd(energy_kcal))

nrow(ukb)
count(ukb[abs(ukb$energy_zscore) > 3, ]) # 268 outliers
ukb <- ukb %>% filter(!abs(energy_zscore) > 3)
nrow(ukb) # 23,746 left

# Check extreme values
ukb[ukb$energy_kcal == min(ukb$energy_kcal),
    c("id", "age", "height", "bmi", "energy_kcal")]
ukb[ukb$energy_kcal == max(ukb$energy_kcal),
    c("id", "age", "height", "bmi", "energy_kcal")]

ggplot(data = ukb) +
  geom_histogram(aes(x = energy_kcal),
                 colour = 1,
                 binwidth = 200) +
  labs(title = "Total energy intake: After trimming",
       y = "Frequency", x = "Total energy intake (kcal)") +
  theme_light()

qqnorm(ukb$energy_kcal, main = "Total energy intake: After trimming")
qqline(ukb$energy_kcal)
# Reasonably normal

## Event-free survival time -----
favstats(ukb$time) # in months


## Graphs by PA level -----

## Average vector magnitude
ggplot(data = ukb) +
  geom_boxplot(aes(y = accel_ave, x = accel_cat)) +
  labs(title = "Boxplot: Accelerometer PA",
       y = "Average vector magnitude (mg)",
       x = "Physical activity level") +
  coord_cartesian(ylim = c(0, 110))

## Age
ggplot(data = ukb) +
  geom_boxplot(aes(y = age, x = accel_cat)) +
  labs(title = "Boxplot: Age",
       y = "Age (years)",
       x = "Physical activity level")

## FEV1
ggplot(data = ukb) +
  geom_boxplot(aes(y = fev1, x = accel_cat)) +
  labs(title = "Boxplot: FEV1",
       y = "FEV1 (L/min)",
       x = "Physical activity level")

## BMI
ggplot(data = ukb) +
  geom_boxplot(aes(y = bmi, x = accel_cat)) +
  labs(title = "Boxplot: BMI",
       y = "BMI (kg/m^2)",
       x = "Physical activity level")

## Total energy intake
ggplot(data = ukb) +
  geom_boxplot(aes(y = energy_kcal, x = accel_cat)) +
  labs(title = "Total energy intake",
       y = "Total energy intake (kcal)",
       x = "Physical activity level")

## Log ASMM
ggplot(data = ukb) +
  geom_boxplot(aes(y = asmm_log, x = accel_cat)) +
  labs(title = "Log of appendicular skeletal mass (kg)",
       y = "Log of ASMM (kg)",
       x = "Physical activity level")

# CATEGORICAL / BINARY -----

## Collapse smoking and alcohol to 2 categories only: Ever vs Never
# lavaan cannot yet handle nominal variables
freq(ukb$smoking)
freq(ukb$alcohol)

# 3 levels
ukb$smoking2 <- ukb$smoking
ukb$alcohol2 <- ukb$alcohol

# 2 levels
ukb <- ukb %>% mutate(
  smoking = fct_collapse(smoking,
                         "Never" = c("Never"),
                         "Ever" = c("Current", "Previous")),
  alcohol = fct_collapse(alcohol,
                         "Never" = c("Never"),
                         "Ever" = c("Current", "Previous")))

## Smoking ----
crosstab(ukb$smoking, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## Alcohol ----
crosstab(ukb$smoking, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## Ethnicity ----
crosstab(ukb$ethnicity, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## Deprivation category -----
ukb$dep_cat <- as.factor(ukb$dep_cat) # Convert to categorical
glimpse(ukb$dep_cat)

ukb %>%
  group_by(dep_cat) %>%
  summarise(MEAN = mean(dep, na.rm=TRUE),
            SD = sd(dep, na.rm=TRUE),
            MEDIAN = median(dep, na.rm=TRUE),
            QUANTILE25 = quantile(dep, c(0.25), na.rm=TRUE),
            QUANTILE75 = quantile(dep, c(0.75), na.rm=TRUE),
            MIN = min(dep, na.rm=TRUE),
            MAX = max(dep, na.rm=TRUE)) %>%
  as.data.frame()
# Townsend index: greater, more deprived
# 1 - least deprived, 10 - most deprived

crosstab(ukb$dep_cat, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## BMI category -----
crosstab(ukb$bmi_cat, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## Sex -----
crosstab(ukb$sex, ukb$accel_cat,
         prop.r = TRUE,
         plot = FALSE)

## Multimorbidity -----
### Charlson comorbidity index ----
charlsonmap <- read.csv("./multimorbidity/charlsoncodes.csv")
names(charlsonmap)[1:2] <- c("ukbcode","ukbdesc")
glimpse(charlsonmap)

ukbdz <- ukb %>%
  select(id, starts_with("illness_"), starts_with("cancer_"))
glimpse(ukbdz)

# Reshape to long format
ukbdz_long <- ukbdz %>%
  pivot_longer(
    cols = c(starts_with("illness_"), starts_with("cancer_")),
    names_to = "position",
    values_to = "ukbcode",
    values_drop_na = TRUE)
print(ukbdz_long)

# Merge ukbdz with charlsonmap
ukbdz_merged <- as_tibble(merge(ukbdz_long, charlsonmap, by = c("ukbcode"),
                                all.x = FALSE, # exclude NA rows
                                all.y = FALSE))
print(ukbdz_merged)
rm(ukbdz, ukbdz_long, charlsonmap)

# Remove duplicate charlsoncode per participant ID
ukbdz_merged[duplicated(ukbdz_merged[, c("id", "charlsoncode")]), ]
ukbdz_merged <- subset(ukbdz_merged,
                       !duplicated(ukbdz_merged[, c("id", "charlsoncode")]))
nrow(ukbdz_merged)

# Create dummy variables for charlsoncodes
ukbdz_merged <- ukbdz_merged %>%
  mutate(
    ci.chf = case_when(charlsoncode == "chf" ~ 1, TRUE ~ 0),
    ci.dementia = case_when(charlsoncode == "dementia" ~ 1, TRUE ~ 0),
    ci.chronicpulm = case_when(charlsoncode == "chronicpulm" ~ 1, TRUE ~ 0),
    ci.rheum = case_when(charlsoncode == "rheum" ~ 1, TRUE ~ 0),
    ci.livermild = case_when(charlsoncode == "livermild" ~ 1, TRUE ~ 0),
    ci.diabetes = case_when(charlsoncode == "diabetes" ~ 1, TRUE ~ 0),
    ci.plegia = case_when(charlsoncode == "plegia" ~ 1, TRUE ~ 0),
    ci.renal = case_when(charlsoncode == "renal" ~ 1, TRUE ~ 0),
    ci.malignancy = case_when(charlsoncode == "malignancy" ~ 1, TRUE ~ 0),
    ci.livermodsev = case_when(charlsoncode == "livermodsev" ~ 1, TRUE ~ 0),
    ci.mets = case_when(charlsoncode == "mets" ~ 1, TRUE ~ 0),
    ci.hivaids = case_when(charlsoncode == "hivaids" ~ 1, TRUE ~ 0))

ukbdz_merged$ukbcode <- NULL
ukbdz_merged$ukbdesc <- NULL
ukbdz_merged <- ukbdz_merged[order(ukbdz_merged$id), ]
print(ukbdz_merged)

# Retain one row per ID
c1 <- aggregate(ukbdz_merged$ci.chf, by = list(id = ukbdz_merged$id), FUN = sum)
c2 <- aggregate(ukbdz_merged$ci.dementia, by = list(id = ukbdz_merged$id), FUN = sum)
c3 <- aggregate(ukbdz_merged$ci.chronicpulm, by = list(id = ukbdz_merged$id), FUN = sum)
c4 <- aggregate(ukbdz_merged$ci.rheum, by = list(id = ukbdz_merged$id), FUN = sum)
c5 <- aggregate(ukbdz_merged$ci.livermild, by = list(id = ukbdz_merged$id), FUN = sum)
c6 <- aggregate(ukbdz_merged$ci.diabetes, by = list(id = ukbdz_merged$id), FUN = sum)
c7 <- aggregate(ukbdz_merged$ci.plegia, by = list(id = ukbdz_merged$id), FUN = sum)
c8 <- aggregate(ukbdz_merged$ci.renal, by = list(id = ukbdz_merged$id), FUN = sum)
c9 <- aggregate(ukbdz_merged$ci.malignancy, by = list(id = ukbdz_merged$id), FUN = sum)
c10 <- aggregate(ukbdz_merged$ci.livermodsev, by = list(id = ukbdz_merged$id), FUN = sum)
c11 <- aggregate(ukbdz_merged$ci.mets, by = list(id = ukbdz_merged$id), FUN = sum)
c12 <- aggregate(ukbdz_merged$ci.hivaids, by = list(id = ukbdz_merged$id), FUN = sum)

names(c1)[2] <- "ci.chf"
names(c2)[2] <- "ci.dementia"
names(c3)[2] <- "ci.chronicpulm"
names(c4)[2] <- "ci.rheum"
names(c5)[2] <- "ci.livermild"
names(c6)[2] <- "ci.diabetes"
names(c7)[2] <- "ci.plegia"
names(c8)[2] <- "ci.renal"
names(c9)[2] <- "ci.malignancy"
names(c10)[2] <- "ci.livermodsev"
names(c11)[2] <- "ci.mets"
names(c12)[2] <- "ci.hivaids"

ukbci <- merge(c1, c2, by = "id")
ukbci <- merge(ukbci, c3, by = "id")
ukbci <- merge(ukbci, c4, by = "id")
ukbci <- merge(ukbci, c5, by = "id")
ukbci <- merge(ukbci, c6, by = "id")
ukbci <- merge(ukbci, c7, by = "id")
ukbci <- merge(ukbci, c8, by = "id")
ukbci <- merge(ukbci, c9, by = "id")
ukbci <- merge(ukbci, c10, by = "id")
ukbci <- merge(ukbci, c11, by = "id")
ukbci <- as_tibble(merge(ukbci, c12, by = "id"))

# Add weights for Charlson score (Quan 2011)
ukbci <- ukbci %>%
  mutate(ci.score =
           ci.chf*2 +ci.dementia*2 + ci.chronicpulm*1 +
           ci.rheum*1 + ci.livermild*2 + ci.diabetes*1 +
           ci.plegia*2 + ci.renal*1 + ci.malignancy*2 +
           ci.livermodsev*4 + ci.mets*6 + ci.hivaids*4)

# Merge with original UKB dataset
ukb <- merge(ukb, ukbci, by = "id", all = TRUE)

### Binary variable for multimorbidity ----
# Replace NA scores with 0
ukb[is.na(ukb$ci.score), c("ci.score")] <- 0
freq(ukb$ci.score)

# Multimorbidity if CCI > 1 (cf. Hanlon 2021)
ukb <- ukb %>%
  mutate(multimorb = case_when(ci.score > 1 ~ 1, TRUE ~ 0))
ukb$multimorb <- as.factor(ukb$multimorb)
levels(ukb$multimorb) <- c("No","Yes")

## Tuberculosis ----
ukb <- ukb %>%
  mutate(tuberculosis = case_when(if_any(.cols = starts_with("illness_"),
                                         .fns = ~. == 1440) ~ 1,
                                  TRUE ~ 0))

# Drop illness and cancer columns
ukb <- ukb %>% select(!c(starts_with("illness_"), starts_with("cancer_")))
glimpse(ukb)


# Summarise by activity level -----
ukb %>%
  group_by(accel_cat) %>%
  summarise(MEAN = mean(accel_ave, na.rm=TRUE),
            SD = sd(accel_ave, na.rm=TRUE),
            MEDIAN = median(accel_ave, na.rm=TRUE),
            QUANTILE25 = quantile(accel_ave, c(0.25), na.rm=TRUE),
            QUANTILE75 = quantile(accel_ave, c(0.75), na.rm=TRUE),
            MIN = min(accel_ave, na.rm=TRUE),
            MAX = max(accel_ave, na.rm=TRUE)) %>%
  as.data.frame()


# BASELINE CHARACTERISTICS ----
## Tabulate all variables ----
tabvars <- c("age", "sex", "ethnicity", "dep_cat",
             "smoking", "alcohol", "energy_kcal",
             "multimorb", "tuberculosis",
             "fev1", "pollution", "spouse",
             "accel_ave", "grip_ave", "bmi_cat", "asmm")
tabnnorm <- c("pollution", "asmm", "accel_ave")
tabcats <- c("bmi_cat", "smoking", "alcohol", "dep_cat", "ethnicity",
             "sex", "multimorb", "tuberculosis", "spouse")

table1 <- CreateTableOne(vars = tabvars,
                         data = ukb,
                         strata = "accel_cat",
                         addOverall = TRUE,
                         factorVars = tabcats,
                         includeNA = TRUE)
table1print <-
  print(table1,
        nonnormal = tabnnorm,
        showAllLevels = TRUE,
        formatOptions = list(big.mark = ","))
# Very few patients with tuberculosis

## Save to a CSV file ----
write.csv(table1print, file = "table1.csv", row.names = TRUE)


# SAVE DATASET WITH FINAL COHORT ----
dim(ukb) # 23,746 rows
saveRDS(ukb, file.path(getwd(), "./UKB_dataset/ukb_finalcohort_analysis.rds"))

