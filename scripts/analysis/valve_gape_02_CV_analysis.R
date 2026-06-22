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



# Model A Figures

# Get model predictions over continuous time (0 to 30 min)
pred <- ggpredict(
  mod_continuous,
  terms = c("time_min [0:30 by=0.5]", "Group", "Status")
)

pred <- as.data.frame(pred)

# Back-transform from log scale to original CV scale
pred <- pred %>%
  mutate(
    predicted = exp(predicted) - 0.01,
    conf.low  = exp(conf.low)  - 0.01,
    conf.high = exp(conf.high) - 0.01
  )

# Figure A1: facet by infection status, color by treatment group
p_con_infection <- ggplot(pred,
                          aes(x = x, y = predicted,
                              color = group, fill = group,
                              group = interaction(group, facet))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),  # ribbon first (behind line)
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +                          # line on top of ribbon
  facet_wrap(~ facet) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  scale_fill_manual(values = c("Control" = "#1f78b4",
                               "CM"      = "#33a02c",
                               "HT"      = "#e31a1c",
                               "HS"      = "#ff7f00"),
                    guide = "none") +
  labs(x = "Time (min)",
       y = "CV of valve gape",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_con_infection

ggsave("outputs/figures/valve_gape_CV_continuous_facet_infection.pdf",
       plot = p_con_infection, width = 6.3, height = 4)

# Figure A2: facet by treatment group, color by infection status
p_con_group <- ggplot(pred,
                      aes(x = x, y = predicted,
                          color = facet, fill = facet,
                          group = interaction(facet, group))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),  # ribbon first
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +                          # line on top
  facet_wrap(~ group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  scale_fill_manual(values = c("uninfected" = "#00BFC4",
                               "infected"   = "#F8766D"),
                    guide = "none") +
  labs(x = "Time (min)",
       y = "CV of valve gape",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_con_group

ggsave("outputs/figures/valve_gape_CV_continuous_facet_group.pdf",
       plot = p_con_group, width = 6.3, height = 4)

# Raw data figure
p_spaghetti_infection <- ggplot(
  total,
  aes(x = time_min,
      y = VC_5min_logCV,
      color = Group,
      group = interaction(Mussel_ID, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "line",
    linewidth = 1.5) +
  facet_wrap(~ Status) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  scale_y_continuous(labels = function(y) round(exp(y) - 0.01, 2)) +
  labs(x = "Time (min)",
       y = "CV of valve gape",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_spaghetti_infection

ggsave("outputs/figures/valve_gape_CV_raw_continuous.pdf",
       plot = p_spaghetti_infection, width = 6.3, height = 4)

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

# Model B Figures

# Figure B1: facet by treatment group, color by infection status
p_cat_group <- ggplot(emm_df,
                      aes(x = time_bin, y = emmean,
                          color = Status, group = Status)) +
  geom_line(linewidth = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2)) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4", "infected" = "#F8766D")) +
  # Restore log-transformed y-axis labels to original CV scale
  scale_y_continuous(
    labels = function(y) round(exp(y) - 0.01, 1)
  ) +
  # Update y-axis label to reflect back-transformed CV values
  labs(y = "CV of valve gape",
       x = "Time bin (5 min)",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_cat_group

ggsave("outputs/figures/valve_gape_CV_categorical_facet_group.pdf",
       plot = p_cat_group, width = 6.3, height = 4)

# Figure B2: facet by infection status, color by treatment group
p_cat_infection <- ggplot(emm_df,
                          aes(x = time_bin, y = emmean,
                              color = Group, group = Group)) +
  geom_line(linewidth = 1.5, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.2)) +
  facet_wrap(~ Status, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  # Restore log-transformed y-axis labels to original CV scale
  scale_y_continuous(
    labels = function(y) round(exp(y) - 0.01, 1)
  ) +
  labs(y = "CV of valve gape",
       x = "Time bin (5 min)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16)
  )

p_cat_infection

ggsave("outputs/figures/valve_gape_CV_categorical_facet_infection.pdf",
       plot = p_cat_infection, width = 6.3, height = 4)

# Raw data figure
p_spaghetti_cat_infection <- ggplot(
  total,
  aes(x = time_bin,
      y = VC_5min_logCV,
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
       y = "CV of valve gape",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_spaghetti_cat_infection

ggsave("outputs/figures/valve_gape_CV_raw_categorical.pdf",
       plot = p_spaghetti_cat_infection, width = 6.3, height = 4)

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