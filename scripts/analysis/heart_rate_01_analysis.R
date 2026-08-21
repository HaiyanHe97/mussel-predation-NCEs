# Heart Rate Analysis
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   Response variable: heart rate (bpm) per 5-min interval
#   Heart rate was estimated using the heartbeatr package. Raw PULSE data were
#   processed with a 30-sec window width and 60-sec window shift. Heart rate
#   estimates within each 5-min window were summarised using the median
#   (pulse_summarise, FUN = median), which is robust to occasional abnormal
#   readings within each window caused by signal noise.
#   Fixed effects: Status (infected/uninfected) * Group (Control/CM/HT/HS) * time
#   Random effects: Rep / Beaker_ID (nested)
#   Two approaches for time:
#     - Model A: time as continuous variable (to assess trends over time)
#     - Model B: time as categorical variable (to assess differences at
#       each time point)
#
# Individual exclusion criteria:
#   23 individuals were excluded prior to statistical analysis
#   (Mussel_5, Mussel_11, Mussel_33, Mussel_45, Mussel_48, Mussel_73,
#   Mussel_88, Mussel_115, Mussel_127, Mussel_143, Mussel_148, Mussel_153,
#   Mussel_166, Mussel_167, Mussel_193, Mussel_194, Mussel_203, Mussel_213,
#   Mussel_224, Mussel_243, Mussel_246, Mussel_283, Mussel_285).
#   Exclusion criterion: highly irregular signal patterns with no discernible
#   rhythmic periodicity. Visual inspection confirmed that the signal pattern
#   of these channels matched that of empty-beaker controls (no mussel),
#   indicating sensor noise rather than genuine cardiac activity.

# 1. Load libraries

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(ggeffects)
library(car)

# 2. Load and format data

hr_total <- read.csv(
  "data/processed/heartrate_merged_final.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
with(unique(hr_total[, c("id", "Group", "Status")]), table(Group, Status))

# Remove flagged individuals
# Criterion: irregular signal pattern matching empty-beaker controls,
# indicating sensor noise rather than genuine cardiac activity
remove_mussels <- c(
  "Mussel_5", "Mussel_11", "Mussel_33", "Mussel_45", "Mussel_48",
  "Mussel_73", "Mussel_88", "Mussel_115", "Mussel_127", "Mussel_143",
  "Mussel_148", "Mussel_153", "Mussel_166", "Mussel_167", "Mussel_193",
  "Mussel_194", "Mussel_203", "Mussel_213", "Mussel_224", "Mussel_243",
  "Mussel_246", "Mussel_283", "Mussel_285"
)

hr_total <- subset(hr_total, !(id %in% remove_mussels))

# Keep time bins 1-5
# Recording duration varied across individuals. In the raw data,
# 24 individuals have records only up to bin 5 (25 min), with no
# missing bpm values and normal heart rate values (mean = 27.5 bpm),
# indicating that their recordings genuinely ended at 25 min rather
# than reflecting signal quality issues. The remaining individuals
# were recorded for longer (up to bin 8). Only bins 1-5 are retained
# to ensure all individuals have complete data across all time points,
# avoiding unbalanced sample sizes in the statistical model.

hr_total$i <- as.numeric(as.character(hr_total$i))
hr_total   <- subset(hr_total, i %in% 1:5)

# Format shared variables
hr_total$Group      <- factor(hr_total$Group,   levels = c("Control", "CM", "HT", "HS"))
hr_total$Status     <- factor(hr_total$Status,  levels = c("uninfected", "infected"))
hr_total$id         <- factor(hr_total$id)
hr_total$Rep        <- factor(hr_total$Rep)
hr_total$Beaker_ID  <- factor(hr_total$Beaker_ID)
hr_total$bpm        <- as.numeric(hr_total$bpm)
hr_total$Length.mm. <- as.numeric(hr_total$Length.mm.)
hr_total$Intensity  <- as.numeric(hr_total$Intensity)

# Check distribution
hist(hr_total$bpm,      main = "Original bpm",        xlab = "bpm")
hist(log(hr_total$bpm), main = "Log-transformed bpm", xlab = "log(bpm)")

# Check sample sizes
table(hr_total$Group, hr_total$Status)
with(unique(hr_total[, c("id", "Group", "Status")]), table(Group, Status))
View(hr_total)

# 3. MODEL A: Time as continuous variable
#    Purpose: assess whether heart rate changes over time (trends)
#    and whether these trends differ by Group and/or infection Status

# Prepare continuous time variable
hr_total$time_min <- as.numeric(as.character(hr_total$i)) * 5 - 2.5

# Fit Model A
mod_hr_continuous <- lmer(
  bpm ~ Status * Group * time_min + (1 | Rep/Beaker_ID),
  data = hr_total
)

summary(mod_hr_continuous)

# Type III ANOVA
anova(mod_hr_continuous)

# Likelihood ratio tests
drop1(mod_hr_continuous, test = "Chisq")

# Model diagnostics
plot(mod_hr_continuous)
qqnorm(resid(mod_hr_continuous)); qqline(resid(mod_hr_continuous))
check_model(mod_hr_continuous)

# Post-hoc: time slopes
emtrends(mod_hr_continuous,
         ~ Group * Status,
         var = "time_min")

emtrends(mod_hr_continuous,
         pairwise ~ Group | Status,
         var = "time_min")

emtrends(mod_hr_continuous,
         pairwise ~ Status | Group,
         var = "time_min")



# 5. MODEL B: Time as categorical variable
#    Purpose: assess heart rate differences at each specific time point
#    and identify which time points show significant group/status differences

# Prepare categorical time variable
hr_total$time_bin <- factor(hr_total$i)

# Fit Model B
mod_hr_categorical <- lmer(
  bpm ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
  data = hr_total
)

summary(mod_hr_categorical)

# Type III ANOVA
anova(mod_hr_categorical)

# Model diagnostics
plot(mod_hr_categorical)
qqnorm(resid(mod_hr_categorical)); qqline(resid(mod_hr_categorical))
check_model(mod_hr_categorical)
simulateResiduals(mod_hr_categorical) |> plot()

# Post-hoc comparisons
# Group differences within each Status x time_bin
emm_hr_group <- emmeans(mod_hr_categorical, ~ Group | Status * time_bin)
pairs(emm_hr_group)

# Status differences within each Group x time_bin
emm_hr_status <- emmeans(mod_hr_categorical, ~ Status | Group * time_bin)
pairs(emm_hr_status)



# 7. Table S0: Shell length and infection intensity as covariates

# Model A versions (time continuous)
mod_hr_length_cont <- lmer(
  bpm ~ Status * Group * time_min + Length.mm. + (1 | Rep/Beaker_ID),
  data = hr_total
)
summary(mod_hr_length_cont)
anova(mod_hr_length_cont)

mod_hr_intensity_cont <- lmer(
  bpm ~ Status * Group * time_min + Intensity + (1 | Rep/Beaker_ID),
  data = hr_total
)
summary(mod_hr_intensity_cont)
anova(mod_hr_intensity_cont)

# Model B versions (time categorical)
mod_hr_length_cat <- lmer(
  bpm ~ Status * Group * time_bin + Length.mm. + (1 | Rep/Beaker_ID),
  data = hr_total
)
summary(mod_hr_length_cat)
anova(mod_hr_length_cat)

mod_hr_intensity_cat <- lmer(
  bpm ~ Status * Group * time_bin + Intensity + (1 | Rep/Beaker_ID),
  data = hr_total
)
summary(mod_hr_intensity_cat)
anova(mod_hr_intensity_cat)
