# CV of Valve Gape Analysis
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   Response variable: coefficient of variation (CV) of valve gape per 5-min interval
#   CV is log-transformed (log(CV + 0.01)) to meet normality assumptions
#   Fixed effects: Status (infected/uninfected) * Group (Control/CM/HT/HS) * time
#   Random effects: Rep / Beaker_ID (nested)
#   Two approaches for time:
#     - Model A: time as continuous variable (to assess trends over time)
#     - Model B: time as categorical variable (to assess differences at
#       each time point)
#   Note: LMM used instead of GLMM (Beta) because response is log-transformed
#   CV, not a bounded proportion. No need to standardize time (no
#   convergence issues).

# 1. Load libraries

library(dplyr)
library(lubridate)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(ggeffects)
library(car)

# 2. Load and format data

total <- read.csv(
  "data/processed/gap_5min_CV.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
View(total)

# Format shared variables
total$Group      <- factor(total$Group,      levels = c("Control", "CM", "HT", "HS"))  # predation risk cue treatment
total$Status     <- factor(total$Status,     levels = c("uninfected", "infected"))      # infection status
total$Mussel_ID  <- factor(total$Mussel_ID)                                             # individual mussel ID
total$Rep        <- factor(total$Rep)                                                   # experimental replicate (block)
total$Beaker_ID  <- factor(total$Beaker_ID)                                            # beaker ID (preparation order)
total$VC_5min_CV    <- as.numeric(total$VC_5min_CV)                                    # response: CV of valve gape
total$Length.mm.    <- as.numeric(gsub(",", ".", total$Length.mm.))                    # shell length (mm)
total$Intensity     <- as.numeric(total$Intensity)                                     # infection intensity (parasite count)

# Log-transform CV to meet normality assumptions
# Adding 0.01 to avoid log(0) for cases where CV = 0
total <- total %>%
  mutate(VC_5min_logCV = log(VC_5min_CV + 0.01))

# Check distribution
hist(total$VC_5min_CV,    main = "Original CV",        xlab = "CV")
hist(total$VC_5min_logCV, main = "Log-transformed CV", xlab = "log(CV + 0.01)")

# Check sample sizes
table(total$Group, total$Status)

# 3. MODEL A: Time as continuous variable
#    Purpose: assess whether CV of valve gape changes over time (trends)
#    and whether these trends differ by Group and/or infection Status

# Prepare continuous time variable
total$time_min <- as.numeric(total$time_mid_min)  # midpoint of each 5-min bin (minutes)
# Note: time is NOT standardized here because LMM does not have the same
# convergence issues as GLMM with Beta family

# Fit Model A
mod_continuous <- lmer(
  VC_5min_logCV ~ Status * Group * time_min + (1 | Rep/Beaker_ID),
  data = total
)

summary(mod_continuous)

# Type III ANOVA
anova(mod_continuous)

# Likelihood ratio tests
drop1(mod_continuous, test = "Chisq")

# Model diagnostics
plot(mod_continuous)                                    # residuals vs fitted
qqnorm(resid(mod_continuous)); qqline(resid(mod_continuous))  # normality check
check_model(mod_continuous)                             # visual diagnostics (performance)

# Post-hoc: compare time slopes across Group * Status combinations
# emtrends compares the slope of time_min for each combination
# Negative slope = CV decreases over time (mussel activity becomes more regular)
emtrends(mod_continuous,
         ~ Group * Status,
         var = "time_min")

# Pairwise slope comparisons between Groups, within each Status level
emtrends(mod_continuous,
         pairwise ~ Group | Status,
         var = "time_min")

# Pairwise slope comparisons between Status levels, within each Group
emtrends(mod_continuous,
         pairwise ~ Status | Group,
         var = "time_min")




# 4. MODEL B: Time as categorical variable
#    Purpose: assess CV differences at each specific time point
#    and identify which time points show significant group/status differences

# Prepare categorical time variable
total$time_bin <- factor(total$time_bin)  # 5-min bins as categories (1-6)

# Fit Model B
mod_categorical <- lmer(
  VC_5min_logCV ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
  data = total
)

summary(mod_categorical)

# Type III ANOVA
anova(mod_categorical)

# Model diagnostics
plot(mod_categorical)
qqnorm(resid(mod_categorical)); qqline(resid(mod_categorical))
check_model(mod_categorical)

# Post-hoc comparisons
# Group differences within each Status x time_bin combination
emm_group <- emmeans(mod_categorical, ~ Group | Status * time_bin)
pairs(emm_group)

# Status differences within each Group x time_bin combination
emm_status <- emmeans(mod_categorical, ~ Status | Group * time_bin)
pairs(emm_status)

# Get model predictions for plotting
emm_pred <- emmeans(mod_categorical, ~ Status | Group * time_bin)
emm_df   <- as.data.frame(emm_pred)
emm_df$time_bin <- factor(emm_df$time_bin)


# 5. Supplementary: Shell length and infection intensity
#    Neither had a significant effect on CV of valve gape

# Model with shell length and infection intensity as covariates
# Included to verify that individual variation in size or parasite load
# does not confound the main treatment effects
model_CV_full <- lmer(
  VC_5min_logCV ~ Group * Status * time_bin + Length.mm. + Intensity +
    (1 | Rep),
  data = total
)
summary(model_CV_full)
anova(model_CV_full)
# Result: neither Length.mm. nor Intensity were significant
# Final models therefore exclude these covariates