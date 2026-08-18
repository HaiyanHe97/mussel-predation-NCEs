# ============================================================
# Figure 4: CV of valve gape
#   (a) time continuous  (facet by infection status, colour by cue)
#   (b) time categorical (5-min bins; facet by status, colour by cue)
#   Figure S7: raw CV trajectories (spaghetti)
#
# Structure: ANALYSIS block fits models + extracts predictions;
#            FIGURE block plots. All fonts bold, size 16.
# Response modelled on log(CV + 0.01); predictions back-transformed
# to the original CV scale for plotting.
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)

# ---- Shared settings ----
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
# ANALYSIS BLOCK
# ============================================================

# ---- Load and format ----
total <- read.csv("data/processed/gap_5min_CV.csv",
                  header = TRUE, sep = ",",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")

total$Group      <- factor(total$Group, levels = c("Control", "CM", "HT", "HS"))
total$Status     <- factor(total$Status, levels = c("uninfected", "infected"))
total$Mussel_ID  <- factor(total$Mussel_ID)
total$Rep        <- factor(total$Rep)
total$Beaker_ID  <- factor(total$Beaker_ID)
total$VC_5min_CV <- as.numeric(total$VC_5min_CV)
total$time_min   <- as.numeric(total$time_mid_min)
total$time_bin   <- factor(total$time_bin)
total <- total %>% mutate(VC_5min_logCV = log(VC_5min_CV + 0.01))

# ---- Model A: time continuous (LMM on log CV) ----
mod_continuous <- lmer(
  VC_5min_logCV ~ Status * Group * time_min + (1 | Rep/Beaker_ID),
  data = total)

pred_a <- ggpredict(mod_continuous,
                    terms = c("time_min [0:30 by=0.5]", "Group", "Status"))
pred_a <- as.data.frame(pred_a)

summary(pred_a$predicted)
summary(total$VC_5min_CV)
# Back-transform predictions from log scale to original CV scale
pred_a <- pred_a %>%
  mutate(predicted = exp(predicted) - 0.01,
         conf.low  = exp(conf.low)  - 0.01,
         conf.high = exp(conf.high) - 0.01)

# Make colour variable identical to Model B (name + levels) for legend merge
pred_a <- pred_a %>% rename(Group = group)
pred_a$Group <- factor(pred_a$Group, levels = c("Control", "CM", "HT", "HS"))

# ---- Model B: time categorical (LMM on log CV) ----
mod_categorical <- lmer(
  VC_5min_logCV ~ Status * Group * time_bin + (1 | Rep/Beaker_ID),
  data = total)

emm_b  <- emmeans(mod_categorical, ~ Group | Status * time_bin)
pred_b <- as.data.frame(emm_b)
pred_b$time_bin <- factor(pred_b$time_bin)

# Back-transform emmeans (and SE band) from log scale to original CV scale
pred_b <- pred_b %>%
  mutate(cv      = exp(emmean) - 0.01,
         cv_low  = exp(emmean - SE) - 0.01,
         cv_high = exp(emmean + SE) - 0.01)

# ============================================================
# FIGURE BLOCK
# ============================================================

# ---- Figure 4a: continuous time ----
p_a <- ggplot(pred_a, aes(x = x, y = predicted,
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
  labs(y = "CV of valve gape", x = "Time (min)",
       color = "Predation risk cue", tag = "a")

# ---- Figure 4b: categorical time ----
#   Back-transformed CV values; y-axis title removed (shares panel a's)
p_b <- ggplot(pred_b, aes(x = time_bin, y = cv,
                          color = Group, group = Group)) +
  geom_line(linewidth = 1, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = cv_low, ymax = cv_high),
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

# ---- Combine Figure 4 ----
fig4 <- p_a + p_b + plot_layout(guides = "collect")
fig4
ggsave("outputs/figures/Figure4_CV_valve_gape.pdf",
       plot = fig4, width = 12, height = 5)

# ---- Figure S7: raw CV trajectories (spaghetti) ----
#   Plotted on log scale, y-axis labels back-transformed to CV
figS7 <- ggplot(
  total,
  aes(x = time_min, y = VC_5min_logCV,
      color = Group,
      group = interaction(Mussel_ID, Status, Group))) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = interaction(Status, Group)),
               fun = mean, geom = "line", linewidth = 1.5) +
  facet_wrap(~ Status) +
  scale_color_manual(values = col_group) +
  scale_y_continuous(labels = function(y) round(exp(y) - 0.01, 2)) +
  base_theme +
  labs(x = "Time (min)", y = "CV of valve gape",
       color = "Predation risk cue")

figS7
ggsave("outputs/figures/FigureS7_CV_raw.pdf",
       plot = figS7, width = 9, height = 5)
