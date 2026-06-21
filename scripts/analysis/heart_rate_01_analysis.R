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

# 4. Figures - Model A

# Get model predictions
pred_hr <- ggpredict(
  mod_hr_continuous,
  terms = c("time_min [0:25 by=0.5]", "Group", "Status")
)

pred_hr <- as.data.frame(pred_hr)

# Model prediction plot: facet by infection status, color by Group
p_hr_con_infection <- ggplot(pred_hr,
                             aes(x = x, y = predicted,
                                 color = group, fill = group,
                                 group = interaction(group, facet))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
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
       y = "Heart rate (bpm)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_con_infection

ggsave("outputs/figures/heart_rate_continuous_facet_infection.pdf",
       plot = p_hr_con_infection, width = 6.3, height = 4)

# Model prediction plot: facet by Group, color by infection status
p_hr_con_group <- ggplot(pred_hr,
                         aes(x = x, y = predicted,
                             color = facet, fill = facet,
                             group = interaction(facet, group))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  scale_fill_manual(values = c("uninfected" = "#00BFC4",
                               "infected"   = "#F8766D"),
                    guide = "none") +
  labs(x = "Time (min)",
       y = "Heart rate (bpm)",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_con_group

ggsave("outputs/figures/heart_rate_continuous_facet_group.pdf",
       plot = p_hr_con_group, width = 6.3, height = 4)

# Figure: Slope dot plot - estimated time slopes by Group x Status
# (Model A post-hoc: emtrends)

# Get slope estimates
slopes_df <- as.data.frame(
  emtrends(mod_hr_continuous,
           ~ Group * Status,
           var = "time_min")
)

# Plot
p_hr_slopes <- ggplot(slopes_df,
                      aes(x = Group,
                          y = time_min.trend,
                          color = Status,
                          shape = Status)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50",
             linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15,
                linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.5)) +
  annotate("text",
           x = 4, y = 0.38,
           label = "*",
           size = 10,
           color = "grey30") +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  scale_shape_manual(values = c("uninfected" = 16,
                                "infected"   = 17)) +
  labs(x     = "Predation risk cue",
       y     = "Heart rate slope",
       color = "Infection status",
       shape = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_slopes

ggsave("outputs/figures/heart_rate_slopes_dotplot.pdf",
       plot = p_hr_slopes, width = 6.3, height = 4)

# Spaghetti plot: facet by infection status, color by Group
p_hr_spaghetti_infection <- ggplot(
  hr_total,
  aes(x = time_min,
      y = bpm,
      color = Group,
      group = interaction(id, Status, Group))) +
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
  labs(x = "Time (min)",
       y = "Heart rate (bpm)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_spaghetti_infection

ggsave("outputs/figures/heart_rate_raw_continuous_facet_infection.pdf",
       plot = p_hr_spaghetti_infection, width = 6.3, height = 4)

# Spaghetti plot: facet by Group, color by infection status
p_hr_spaghetti_group <- ggplot(
  hr_total,
  aes(x = time_min,
      y = bpm,
      color = Status,
      group = interaction(id, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(
    aes(group = interaction(Status, Group)),
    fun = mean,
    geom = "line",
    linewidth = 1.5) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  labs(x = "Time (min)",
       y = "Heart rate (bpm)",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_spaghetti_group

ggsave("outputs/figures/heart_rate_raw_continuous_facet_group.pdf",
       plot = p_hr_spaghetti_group, width = 6.3, height = 4)

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

# 6. Figures - Model B

# Get model predictions
emm_hr_pred <- emmeans(mod_hr_categorical, ~ Group | Status * time_bin)
emm_hr_df   <- as.data.frame(emm_hr_pred)
emm_hr_df$time_bin <- factor(emm_hr_df$time_bin)

# Model prediction plot: facet by Group, color by infection status
p_hr_cat_group <- ggplot(emm_hr_df,
                         aes(x = time_bin, y = emmean,
                             color = Status, group = Status)) +
  geom_line(linewidth = 1,
            position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15,
                position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.2)) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  labs(y = "Heart rate (bpm)",
       x = "Time bin (5 min)",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_cat_group

ggsave("outputs/figures/heart_rate_categorical_facet_group.pdf",
       plot = p_hr_cat_group, width = 6.3, height = 4)

# Model prediction plot: facet by infection status, color by Group
p_hr_cat_infection <- ggplot(emm_hr_df,
                             aes(x = time_bin, y = emmean,
                                 color = Group, group = Group)) +
  geom_line(linewidth = 1.5,
            position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15,
                position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.2)) +
  facet_wrap(~ Status, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#1f78b4",
                                "CM"      = "#33a02c",
                                "HT"      = "#e31a1c",
                                "HS"      = "#ff7f00")) +
  labs(y = "Heart rate (bpm)",
       x = "Time bin (5 min)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_cat_infection

ggsave("outputs/figures/heart_rate_categorical_facet_infection.pdf",
       plot = p_hr_cat_infection, width = 6.3, height = 4)

# Spaghetti plot: facet by infection status, color by Group
p_hr_spaghetti_cat_infection <- ggplot(
  hr_total,
  aes(x = time_bin,
      y = bpm,
      color = Group,
      group = interaction(id, Status, Group))) +
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
  labs(x = "Time bin (5 min)",
       y = "Heart rate (bpm)",
       color = "Predation risk cue") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_spaghetti_cat_infection

ggsave("outputs/figures/heart_rate_raw_categorical_facet_infection.pdf",
       plot = p_hr_spaghetti_cat_infection, width = 6.3, height = 4)

# Spaghetti plot: facet by Group, color by infection status
p_hr_spaghetti_cat_group <- ggplot(
  hr_total,
  aes(x = time_bin,
      y = bpm,
      color = Status,
      group = interaction(id, Status, Group))) +
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
  facet_wrap(~ Group) +
  theme_bw() +
  scale_color_manual(values = c("uninfected" = "#00BFC4",
                                "infected"   = "#F8766D")) +
  labs(x = "Time bin (5 min)",
       y = "Heart rate (bpm)",
       color = "Infection status") +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 16))

p_hr_spaghetti_cat_group

ggsave("outputs/figures/heart_rate_raw_categorical_facet_group.pdf",
       plot = p_hr_spaghetti_cat_group, width = 6.3, height = 4)

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
