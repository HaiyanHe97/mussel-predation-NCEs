# =============================================================
# Supporting analyses
# Mussel predation-risk experiment (infection x predator cue)
#
# What this script does:
#   (1) Summarises effect magnitudes - medians per treatment for each
#       metric, so results can report how much higher/lower a treatment
#       is, not only which differences are significant (Sections 1.2, 2.1).
#   (2) Checks heart-rate data quality - how many 0 / near-0 bpm
#       readings occur, which may be artefacts (e.g. missed beats,
#       mostly in the first time bin) (Section 1.1).
#   (3) Fits the heart-rate mixed models (Section 1.3).
#
# Inputs (prepared upstream, not created here):
#   hr_total   : per-individual heart-rate records.
#                Columns used: id, Status, Group, Rep, Beaker_ID,
#                              i (time bin, 1 = first bin), bpm
#   move_total : per-individual movement / aggregation summaries.
#                Columns used: Status, Group, mean_A, total_time_agg,
#                              Amax, Start_Time, mean_gross_mm,
#                              mean_net_mm, mean_CI, Byssus
#
#   Status : uninfected / infected
#   Group  : Control / CM (native crab) / HT, HS (invasive crabs)
# =============================================================

## ---- Packages ----
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)


# =============================================================
# 1. HEART RATE
# =============================================================

## ---- 1.1 Data-quality check: 0 / near-0 bpm readings -------
##      flag possibly artefactual readings

# Overall range of bpm
summary(hr_total$bpm)

# Count of zero / near-zero readings
hr_total %>% summarise(
  n             = n(),
  n_zero        = sum(bpm == 0,  na.rm = TRUE),
  n_near_zero   = sum(bpm <= 5,  na.rm = TRUE),
  pct_near_zero = mean(bpm <= 5, na.rm = TRUE) * 100
)

# Are low values concentrated in the first time bin?
hr_total %>% group_by(i) %>%
  summarise(mean_bpm      = mean(bpm, na.rm = TRUE),
            min_bpm       = min(bpm,  na.rm = TRUE),
            pct_near_zero = mean(bpm <= 5, na.rm = TRUE) * 100,
            .groups = "drop") %>%
  arrange(i)

# Per-individual: lowest value, bin-1 starting value, any zeros
hr_total %>% group_by(id) %>%
  summarise(min_bpm   = min(bpm, na.rm = TRUE),
            start_bpm = bpm[i == 1][1],
            ever_zero = any(bpm == 0, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(min_bpm) %>% head(20)


## ---- 1.2 Descriptive summaries (effect magnitudes) ---------

# Mean bpm by infection status and cue
hr_total %>%
  group_by(Status, Group) %>%
  summarise(mean_bpm = round(mean(bpm, na.rm = TRUE), 1),
            sd       = round(sd(bpm,   na.rm = TRUE), 1),
            n        = n_distinct(id),
            .groups  = "drop") %>%
  arrange(Status, Group)

# Mean bpm by status, cue and time bin
hr_total %>%
  group_by(Status, Group, i) %>%
  summarise(mean_bpm = round(mean(bpm, na.rm = TRUE), 1),
            se       = round(sd(bpm, na.rm = TRUE) / sqrt(n()), 2),
            .groups  = "drop") %>%
  arrange(Status, Group, i) %>%
  print(n = 40)


## ---- 1.3 Models (CORE) -------------------------------------
##      First time bin removed (see QC in 1.1)

hr_noT0 <- hr_total %>% filter(i != 1)
hr_noT0$time_min <- as.numeric(as.character(hr_noT0$i)) * 5 - 2.5
hr_noT0$time_bin <- droplevels(factor(hr_noT0$i))

# Model A - continuous time
# NOTE: in the original snippets this model was referred to both as
# `modA_noT0` and `mod_hr_continuous`. Assumed to be the same model
# and unified here as `mod_hr_cont`. Confirm before use.
mod_hr_cont <- lmer(bpm ~ Status * Group * time_min + (1 | Rep/Beaker_ID),
                    data = hr_noT0)
anova(mod_hr_cont)
emtrends(mod_hr_cont, pairwise ~ Group  | Status, var = "time_min")  # slope diffs among cues
emtrends(mod_hr_cont, pairwise ~ Status | Group,  var = "time_min")  # slope diffs by infection

# Model B - categorical time
mod_hr_cat <- lmer(bpm ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
                   data = hr_noT0)
anova(mod_hr_cat)
pairs(emmeans(mod_hr_cat, ~ Group | Status * time_bin))  # cue diffs at each time bin


## ---- 1.4 Change over the recording interval ----------------
##      start (2.5 min) vs end (22.5 min); magnitude of bpm change

emmeans(mod_hr_cont, ~ Group | Status,
        at = list(time_min = c(2.5, 22.5)))
emmeans(mod_hr_cont, ~ time_min | Group * Status,
        at = list(time_min = c(2.5, 22.5)))

# Uninfected vs infected at the end point, within each cue
emm_end <- emmeans(mod_hr_cont, ~ Group * Status, at = list(time_min = 22.5))
pairs(emm_end, by = "Group")


# =============================================================
# 2. MOVEMENT & AGGREGATION
# =============================================================

## ---- 2.1 Descriptive summaries (effect magnitudes) ---------

# Infection main effect (pooled across cues)
move_total %>% group_by(Status) %>%
  summarise(med_meanA  = round(median(mean_A,         na.rm = TRUE), 1),
            med_totagg = round(median(total_time_agg, na.rm = TRUE), 1),
            .groups = "drop")

# Cue main effect (pooled across infection status)
move_total %>% group_by(Group) %>%
  summarise(med_meanA  = round(median(mean_A,         na.rm = TRUE), 1),
            med_totagg = round(median(total_time_agg, na.rm = TRUE), 1),
            .groups = "drop")

# Median of each metric by cue x status
aggregate(Start_Time    ~ Group + Status, data = move_total, FUN = median)
aggregate(mean_gross_mm ~ Group + Status, data = move_total, FUN = median)
aggregate(mean_net_mm   ~ Group + Status, data = move_total, FUN = median)
aggregate(Byssus        ~ Group + Status, data = move_total, FUN = median)
aggregate(Start_Time    ~ Status,         data = move_total, FUN = median)
aggregate(Amax          ~ Group,          data = move_total, FUN = median)
aggregate(mean_CI       ~ Group,          data = move_total, FUN = median)


## ---- 2.2 Exploratory check: HT, uninfected distance --------
##      small n per subgroup; used to inspect spread / outliers

move_total %>% filter(Group == "HT", Status == "uninfected") %>%
  summarise(n      = n(),
            min    = round(min(mean_gross_mm,            na.rm = TRUE), 1),
            Q1     = round(quantile(mean_gross_mm, 0.25, na.rm = TRUE), 1),
            median = round(median(mean_gross_mm,         na.rm = TRUE), 1),
            Q3     = round(quantile(mean_gross_mm, 0.75, na.rm = TRUE), 1),
            max    = round(max(mean_gross_mm,            na.rm = TRUE), 1))

# The 6 raw values behind that subgroup
move_total %>% filter(Group == "HT", Status == "uninfected") %>%
  pull(mean_gross_mm)
