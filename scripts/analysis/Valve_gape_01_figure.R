# ============================================================
# Valve Gape Figures
# Project: Effects of parasitic infection on predation risk
#          responses in blue mussels (Mytilus edulis)
#
#   Figure 3:  (a) valve gape, time continuous  (facet by infection
#              status, colour by predation risk cue)
#              (b) valve gape, time categorical (5-min bins; facet by
#              infection status, colour by predation risk cue)
#   Figure S6: raw valve gape trajectories (spaghetti plot)
#
# Structure: ANALYSIS block fits models and extracts predictions;
#            FIGURE block uses only those predictions to plot.
# All fonts bold, size 16.
# ============================================================

# ------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(emmeans)
library(ggeffects)

# ------------------------------------------------------------
# SHARED SETTINGS (colours + theme, reused by all panels)
# ------------------------------------------------------------
col_group <- c("Control" = "#1f78b4", "CM" = "#33a02c",
               "HT" = "#e31a1c", "HS" = "#ff7f00")

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
# ANALYSIS BLOCK  (models + predictions only, no plotting)
# ============================================================

# ---- Load and format data ----
total <- read.csv("data/processed/gap_5min_activity.csv",
                  header = TRUE, sep = ",",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")

total$Group      <- factor(total$Group, levels = c("Control", "CM", "HT", "HS"))
total$Status     <- factor(total$Status, levels = c("uninfected", "infected"))
total$Mussel_ID  <- factor(total$Mussel_ID)
total$Rep        <- factor(total$Rep)
total$Beaker_ID  <- factor(total$Beaker_ID)
total$VC_5min_avg <- as.numeric(total$VC_5min_avg)
total$time_min   <- as.numeric(total$time_mid_min)
total$time_min_z <- scale(total$time_min)
total$time_bin   <- factor(total$time_bin)

time_mean <- mean(total$time_min)
time_sd   <- sd(total$time_min)

# ---- Model A: time continuous (Beta GLMM) ----
mod_continuous <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_min_z + (1 | Rep/Beaker_ID),
  family = beta_family(), data = total)

# Predictions at the 6 real time points (5,10,...,30 min)
pred_a <- ggpredict(
  mod_continuous,
  terms = c("time_min_z [-1.46,-0.88,-0.29,0.29,0.88,1.46]", "Group", "Status"))
pred_a <- pred_a %>% arrange(group, facet, x)

# Map z-scores back to minutes for a readable x-axis
time_map <- c("-1.46" = 5, "-0.88" = 10, "-0.29" = 15,
              "0.29" = 20, "0.88" = 25, "1.46" = 30)
pred_a$time_min <- time_map[as.character(round(pred_a$x, 2))]
pred_a$time_min <- factor(pred_a$time_min, levels = c(5, 10, 15, 20, 25, 30))

# Make the colour variable identical to Model B (name + levels) so the
# legend can be merged: rename 'group' -> 'Group', set the same levels
pred_a <- pred_a %>% rename(Group = group)
pred_a$Group <- factor(pred_a$Group, levels = c("Control", "CM", "HT", "HS"))

# ---- Model B: time categorical (Beta GLMM) ----
mod_categorical <- glmmTMB(
  VC_5min_avg ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
  family = beta_family(link = "logit"), data = total)

# Predicted means (response scale) per Group x Status x time_bin
emm_b  <- emmeans(mod_categorical, ~ Group | Status * time_bin,
                  type = "response")
pred_b <- as.data.frame(emm_b)

# ============================================================
# FIGURE BLOCK  (uses pred_a / pred_b / total only)
# ============================================================

# ---- Figure 3a: continuous time ----
#   Keeps the y-axis title; fill legend hidden so only one colour
#   legend remains (identical to panel b -> legends merge)
p_a <- ggplot(pred_a, aes(x = time_min, y = predicted,
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
  labs(y = "Valve gape (fraction open)", x = "Time (min)",
       color = "Predation risk cue", tag = "a")

# ---- Figure 3b: categorical time ----
#   Y-axis title removed (shares panel a's); same colour legend as a
p_b <- ggplot(pred_b, aes(x = time_bin, y = response,
                          color = Group, group = Group)) +
  geom_line(linewidth = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.15, position = position_dodge(width = 0.2)) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.05,
                                             dodge.width = 0.2)) +
  facet_wrap(~ Status, ncol = 1) +
  scale_color_manual(values = col_group) +
  base_theme +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  labs(y = NULL, x = "Time bin (5 min)",
       color = "Predation risk cue", tag = "b")

# ---- Combine Figure 3: a | b, single shared legend on the right ----
fig3 <- p_a + p_b + plot_layout(guides = "collect")
fig3
ggsave("outputs/figures/Figure3_valve_gape.pdf",
       plot = fig3, width = 12, height = 5)

# ---- Figure S6: raw valve gape trajectories (spaghetti) ----
figS6 <- ggplot(
  total,
  aes(x = time_min, y = VC_5min_avg,
      color = Group,
      group = interaction(Mussel_ID, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = interaction(Status, Group)),
               fun = mean, geom = "line", linewidth = 1.5) +
  facet_wrap(~ Status) +
  scale_color_manual(values = col_group) +
  base_theme +
  labs(x = "Time (min)", y = "Valve gape (fraction open)",
       color = "Predation risk cue")

figS6
ggsave("outputs/figures/FigureS6_valve_gape_raw.pdf",
       plot = figS6, width = 9, height = 5)
