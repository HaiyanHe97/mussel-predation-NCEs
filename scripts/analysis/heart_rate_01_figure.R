# ============================================================
# Figure S5: Heart rate (appendix)
#   (a) heart rate over time (continuous; facet by infection status,
#       colour by predation risk cue) + 95% CI
#   (b) estimated temporal slope per cue, by infection status,
#       with 95% CI (dot plot; asterisk marks HS status difference)
#   Figure S8: raw heart rate trajectories (spaghetti)
#
# ANALYSIS block fits models + extracts predictions/slopes;
# FIGURE block plots. All fonts bold, size 16.
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)

# ---- Shared settings ----
col_group  <- c("Control" = "#1f78b4", "CM" = "#33a02c",
                "HT" = "#e31a1c", "HS" = "#ff7f00")
col_status <- c("uninfected" = "#00BFC4", "infected" = "#F8766D")

base_theme <- theme_bw() +
  theme(
    strip.text   = element_text(face = "bold", size = 16, color = "black"),
    axis.text    = element_text(face = "bold", size = 16, color = "black"),
    axis.title   = element_text(face = "bold", size = 16, color = "black"),
    legend.title = element_text(face = "bold", size = 16, color = "black"),
    legend.text  = element_text(face = "bold", size = 16, color = "black"),
    plot.tag     = element_text(face = "bold", size = 16, color = "black"),
    legend.key.size = unit(0.5, "cm")
  )

# ============================================================
# ANALYSIS BLOCK
# ============================================================

# ---- Load, exclude flagged individuals, keep bins 1-5 ----
hr_total <- read.csv("data/processed/heartrate_merged_final.csv",
                     header = TRUE, sep = ",",
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")

remove_mussels <- c(
  "Mussel_5","Mussel_11","Mussel_33","Mussel_45","Mussel_48",
  "Mussel_73","Mussel_88","Mussel_115","Mussel_127","Mussel_143",
  "Mussel_148","Mussel_153","Mussel_166","Mussel_167","Mussel_193",
  "Mussel_194","Mussel_203","Mussel_213","Mussel_224","Mussel_243",
  "Mussel_246","Mussel_283","Mussel_285")
hr_total <- subset(hr_total, !(id %in% remove_mussels))

hr_total$i <- as.numeric(as.character(hr_total$i))
hr_total   <- subset(hr_total, i %in% 1:5)

hr_total$Group     <- factor(hr_total$Group, levels = c("Control","CM","HT","HS"))
hr_total$Status    <- factor(hr_total$Status, levels = c("uninfected","infected"))
hr_total$id        <- factor(hr_total$id)
hr_total$Rep       <- factor(hr_total$Rep)
hr_total$Beaker_ID <- factor(hr_total$Beaker_ID)
hr_total$bpm       <- as.numeric(hr_total$bpm)
hr_total$time_min  <- as.numeric(as.character(hr_total$i)) * 5 - 2.5

# ---- Model A: time continuous (LMM) ----
mod_hr <- lmer(bpm ~ Status * Group * time_min + (1 | Rep/Beaker_ID),
               data = hr_total)

# Predictions over time
pred_a <- as.data.frame(
  ggpredict(mod_hr, terms = c("time_min [0:25 by=0.5]", "Group", "Status")))
# Match colour variable to col_group (name + levels)
pred_a <- pred_a %>% rename(Group = group)
pred_a$Group <- factor(pred_a$Group, levels = c("Control","CM","HT","HS"))

# Slopes per Group x Status (for panel b)
slopes_df <- as.data.frame(
  emtrends(mod_hr, ~ Group * Status, var = "time_min"))

# ============================================================
# FIGURE BLOCK
# ============================================================

# ---- Figure S5a: heart rate over time ----
s5_a <- ggplot(pred_a, aes(x = x, y = predicted,
                           color = Group,
                           group = interaction(Group, facet))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Group),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ facet) +
  scale_color_manual(values = col_group) +
  scale_fill_manual(values = col_group) +
  guides(fill = "none") +
  base_theme +
  labs(y = "Heart rate (bpm)", x = "Time (min)",
       color = "Predation risk cue", tag = "a")

# ---- Figure S5b: temporal slope per cue, by infection status ----
s5_b <- ggplot(slopes_df, aes(x = Group, y = time_min.trend,
                              color = Status, shape = Status)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  annotate("text", x = 4, y = 0.38, label = "*", size = 10, color = "grey30") +
  scale_color_manual(values = col_status) +
  scale_shape_manual(values = c("uninfected" = 16, "infected" = 17)) +
  base_theme +
  labs(x = "Predation risk cue", y = "Heart rate slope (bpm/min)",
       color = "Infection status", shape = "Infection status", tag = "b")

# ---- Combine: a | b (two different legends kept, collected right) ----
figS5 <- s5_a + s5_b + plot_layout(guides = "collect")
figS5
ggsave("outputs/figures/FigureS5_heart_rate.pdf",
       plot = figS5, width = 13, height = 5)

# ---- Figure S8: raw heart rate trajectories (spaghetti) ----
figS8 <- ggplot(hr_total,
                aes(x = time_min, y = bpm, color = Group,
                    group = interaction(id, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = interaction(Status, Group)),
               fun = mean, geom = "line", linewidth = 1.5) +
  facet_wrap(~ Status) +
  scale_color_manual(values = col_group) +
  coord_cartesian(ylim = c(0, 50)) +
  base_theme +
  labs(x = "Time (min)", y = "Heart rate (bpm)",
       color = "Predation risk cue")

figS8
ggsave("outputs/figures/FigureS8_heart_rate_raw.pdf",
       plot = figS8, width = 9, height = 5)
