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

# Model A Figures

# Back-transform the 6 actual time points to z-scores for ggpredict
# z = (x - mean) / sd
# time points: 5, 10, 15, 20, 25, 30 min
z_vals <- round((c(5, 10, 15, 20, 25, 30) - time_mean) / time_sd, 2)
# z_vals = -1.46, -0.88, -0.29, 0.29, 0.88, 1.46

# Get model predictions at actual time points
pred <- ggpredict(
  mod_continuous,
  terms = c("time_min_z [-1.46,-0.88,-0.29,0.29,0.88,1.46]", "Group", "Status")
)

pred <- pred %>% arrange(group, facet, x)

# Map z-scores back to minutes for readable x-axis
time_map <- c(
  "-1.46" = 5,
  "-0.88" = 10,
  "-0.29" = 15,
  "0.29"  = 20,
  "0.88"  = 25,
  "1.46"  = 30
)

pred$time_min <- time_map[as.character(round(pred$x, 2))]
pred$time_min <- factor(pred$time_min, levels = c(5, 10, 15, 20, 25, 30))

# Figure A1: facet by infection status, color by treatment group
p_con_infection <- ggplot(pred, aes(x = time_min, y = predicted,
                                    color = group,
                                    group = interaction(group, facet))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),  # ribbon first (behind line)
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +                                        # line on top of ribbon
  facet_wrap(~ facet) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  scale_fill_manual(values = c("Control" = "#1f78b4",
                               "CM"      = "#33a02c",
                               "HT"      = "#e31a1c",
                               "HS"      = "#ff7f00")) +
  labs(y = "Valve gape (fraction open)",
       x = "Time (min)",
       color = "Predation risk cue",
       fill = "Predation risk cue")    +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_con_infection

ggsave("outputs/figures/valve_gape_continuous_facet_infection.pdf",
       plot = p_con_infection, width = 6.3, height = 4)

# Figure A2: facet by treatment group, color by infection status
p_con_group <- ggplot(pred, aes(x = time_min, y = predicted,
                                color = facet,
                                group = interaction(facet, group))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet),  # ribbon first
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +                                        # line on top
  facet_wrap(~ group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  scale_fill_manual(values = c("uninfected" = "#00BFC4",
                               "infected"   = "#F8766D")) +
  labs(y = "Valve gape (fraction open)",
       x = "Time (min)",
       color = "Infection status",
       fill = "Infection status")   +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_con_group

ggsave("outputs/figures/valve_gape_continuous_facet_group.pdf",
       plot = p_con_group, width = 6.3, height = 4)

# Raw data figure
p_spaghetti_con_infection <- ggplot(
  total,
  aes(x = time_min,
      y = VC_5min_avg,
      color = Group,
      group = interaction(Mussel_ID, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "line",
    linewidth = 1.5) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "line",
    size = 1.5) +
  facet_wrap(~ Status) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  scale_y_continuous(labels = function(y) round(exp(y) - 0.01, 0)) +
  labs(x = "Time (min)",
       y = "Valve gape (fraction open)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_spaghetti_con_infection

ggsave("outputs/figures/valve_gape_raw_continuous.pdf",
       plot = p_spaghetti_con_infection, width = 6.3, height = 4)

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

# Model B Figures

# Figure B1: facet by treatment group, color by infection status
p_cat_group <- ggplot(emm_df, aes(x = time_bin, y = response,
                                  color = Status, group = Status)) +
  geom_line(linewidth = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.15, position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2)) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4", "infected" = "#F8766D")) +
  labs(y = "Valve gape (fraction open)", x = "Time bin (5 min)", color = " Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_cat_group

ggsave("outputs/figures/valve_gape_categorical_facet_group.pdf",
       plot = p_cat_group, width = 6.3, height = 4)

# Figure B2: facet by infection status, color by treatment group
p_cat_infection <- ggplot(emm_df,
                          aes(x = time_bin, y = response,
                              color = Group, group = Group)) +
  geom_line(linewidth = 1.5, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.15, position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2)) +
  facet_wrap(~ Status, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  labs(y = "Valve gape (fraction open)", x = "Time bin (5 min)", color = "Predation risk cue ") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_cat_infection

ggsave("outputs/figures/valve_gape_categorical_facet_infection.pdf",
       plot = p_cat_infection, width = 6.3, height = 4)

# Raw data figure
p_spaghetti_cat_infection <- ggplot(
  total,
  aes(x = time_bin,
      y = VC_5min_avg,
      color = Group,
      group = interaction(Mussel_ID, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "line",
    linewidth = 1.5) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "point",
    size = 2.5) +
  facet_wrap(~ Status) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  scale_y_continuous(labels = function(y) round(exp(y) - 0.01, 2)) +
  labs(x = "Time bin (5 min)",
       y = "Valve gape (fraction open)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_spaghetti_cat_infection

ggsave("outputs/figures/valve_gape_raw_categorical.pdf",
       plot = p_spaghetti_cat_infection, width = 6.3, height = 4)

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
