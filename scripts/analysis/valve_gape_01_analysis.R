# Valve Gape Analysis
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   Response variable: mean fraction openness per 5-min interval (VC_5min_avg)
#   Fixed effects: Status (infected/uninfected) * Group (Control/CM/HT/HS) * time
#   Random effects: Rep / Beaker_ID (nested)
#   Rep: accounts for variation between experimental blocks (batches)
#   Beaker_ID nested within Rep: accounts for preparation order effects -
#   within each block, mussels were prepared in a fixed order (sensor
#   attachment, drying), so earlier-prepared individuals waited longer
#   in beakers before cue water was added, potentially introducing
#   handling-related stress that could influence behavioural responses
#   Two approaches for time:
#     - Model A: time as continuous variable (to assess trends over time)
#     - Model B: time as categorical variable (to assess differences at
#       each time point)

# 1. Load libraries

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(gghalves)
library(mgcv)
library(ggeffects)
library(glmmTMB)
library(car)
library(DHARMa)
library(performance)

# 2. Load and format data

total <- read.csv(
  "data/processed/gap_5min_activity.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
View(total)

# Format shared variables (used in both models)
total$Group      <- factor(total$Group,      levels = c("Control", "CM", "HT", "HS"))  # predation risk cue treatment
total$Status     <- factor(total$Status,     levels = c("uninfected", "infected"))      # infection status
total$Mussel_ID  <- factor(total$Mussel_ID)                                             # individual mussel ID
total$Rep        <- factor(total$Rep)                                                   # experimental replicate (block)
total$Beaker_ID  <- factor(total$Beaker_ID)                                            # beaker ID (mussel placement order)
total$VC_5min_avg   <- as.numeric(total$VC_5min_avg)                                   # response: mean fraction openness (0-1)
total$Length.mm.    <- as.numeric(gsub(",", ".", total$Length.mm.))                    # shell length (mm)
total$Intensity     <- as.numeric(total$Intensity)                                     # infection intensity (parasite count)

# Check distribution of response variable
hist(total$VC_5min_avg, main = "Distribution of valve gape", xlab = "Fraction open")
# Check for exact 0s and 1s (not allowed in Beta distribution)
sum(total$VC_5min_avg == 0, na.rm = TRUE)
sum(total$VC_5min_avg == 1, na.rm = TRUE)

# Check sample sizes
table(total$Group, total$Status)
total %>%
  distinct(Mussel_ID, Group, Status, Intensity) %>%
  filter(Group == "CM") %>%
  arrange(Status, Intensity)

# 3. MODEL A: Time as continuous variable
#    Purpose: assess whether valve gape changes over time (trends)
#    and whether these trends differ by Group and/or infection Status

# Prepare continuous time variable
total$time_min <- as.numeric(total$time_mid_min)  # midpoint of each 5-min bin (minutes)

# Standardize time to improve model convergence in glmmTMB
# (raw time_min causes NA/NaN warnings during optimization)
# Results are identical to non-standardized model (confirmed by AIC)
total$time_min_z <- scale(total$time_min)

# Store mean and SD for back-transforming z-scores to minutes in figures
time_mean <- mean(total$time_min)  # 17.5 min
time_sd   <- sd(total$time_min)    # 8.54 min

# Fit Model A
mod_continuous <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_min_z + (1 | Rep/Beaker_ID),
  family = beta_family(),
  data = total
)

summary(mod_continuous)

# Type III Wald chi-square tests
Anova(mod_continuous, type = 3)

# Likelihood ratio tests
drop1(mod_continuous, test = "Chisq")

# Model diagnostics
check_model(mod_continuous)                       # visual diagnostics (performance)
res_con <- simulateResiduals(mod_continuous)      # simulation-based residuals (DHARMa)
plot(res_con, asFactor = TRUE)
plot(res_con, rank = TRUE)

# Post-hoc: compare time slopes across Group * Status combinations
# emtrends compares the slope of time_min_z for each combination
# Positive slope = valve gape increases over time; negative = decreases
emtrends(mod_continuous,
         ~ Group * Status,
         var = "time_min_z")

# Pairwise slope comparisons between Groups, within each Status level
emtrends(mod_continuous,
         pairwise ~ Group | Status,
         var = "time_min_z")

# Pairwise slope comparisons between Status levels, within each Group
emtrends(mod_continuous,
         pairwise ~ Status | Group,
         var = "time_min_z")



# 4. MODEL B: Time as categorical variable
#    Purpose: assess valve gape differences at each specific time point
#    and identify which time points show significant group/status differences

# Prepare categorical time variable
total$time_bin <- factor(total$time_bin)  # 5-min bins as categories (1-6)

# Fit Model B
mod_categorical <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
  family = beta_family(link = "logit"),
  data = total
)

summary(mod_categorical)

# Type III Wald chi-square tests
Anova(mod_categorical, type = 3)

# Likelihood ratio tests
drop1(mod_categorical, test = "Chisq")

# Model diagnostics
check_model(mod_categorical)
res_cat <- simulateResiduals(mod_categorical)
plot(res_cat, asFactor = TRUE)
plot(res_cat, rank = TRUE)

# Post-hoc comparisons
# Status differences within each Group x time_bin combination
emm_status <- emmeans(mod_categorical, ~ Status | Group * time_bin)
pairs(emm_status)

# Group differences within each Status x time_bin combination
emm_group <- emmeans(mod_categorical, ~ Group | Status * time_bin)
pairs(emm_group)

# Get model predictions for plotting
emm_pred <- emmeans(mod_categorical, ~ Status | Group * time_bin, type = "response")
emm_df   <- as.data.frame(emm_pred)


# 5. Supplementary: Shell length and infection intensity as covariates
#    Valve gape - Table S0

# Model A (time continuous, standardised)

mod_vg_length_cont <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_min_z + Length.mm. + (1 | Rep/Beaker_ID),
  family = beta_family(),
  data = total
)
summary(mod_vg_length_cont)
Anova(mod_vg_length_cont, type = 3)

mod_vg_intensity_cont <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_min_z + Intensity + (1 | Rep/Beaker_ID),
  family = beta_family(),
  data = total
)
summary(mod_vg_intensity_cont)
Anova(mod_vg_intensity_cont, type = 3)

# Model B (time categorical)

mod_vg_length_cat <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_bin + Length.mm. + (1 | Rep/Beaker_ID),
  family = beta_family(link = "logit"),
  data = total
)
summary(mod_vg_length_cat)
Anova(mod_vg_length_cat, type = 3)

mod_vg_intensity_cat <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_bin + Intensity + (1 | Rep/Beaker_ID),
  family = beta_family(link = "logit"),
  data = total
)
summary(mod_vg_intensity_cat)
Anova(mod_vg_intensity_cat, type = 3)
