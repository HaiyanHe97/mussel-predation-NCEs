# # ============================================================
# Aggregation & Movement Figures
# Project: Effects of parasitic infection on predation risk
#          responses in blue mussels (Mytilus edulis)
#
# Main text:
#   Figure 1: (a) mean aggregation proportion, (b) total aggregation time
#   Figure 2: (a) time to first aggregation, (b) mean total distance
# Appendix (ordered by metric number in Table 1):
#   Figure S1: maximum aggregation proportion
#   Figure S2: (a) mean net displacement, (b) mean confinement index
#   Figure S3: byssus thread count
#
# Shared layout: paired panels share the x-axis (Group); a single
# legend is collected to the right; all fonts bold, size 13.
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ---- Data ----
move_total <- read.csv(
  "data/processed/merged_movement.csv",
  header = TRUE, sep = ",", stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
move_total$Group  <- factor(move_total$Group,
                            levels = c("Control", "CM", "HT", "HS"))
move_total$Status <- factor(move_total$Status,
                            levels = c("uninfected", "infected"))

# ---- Colours ----
col_status <- c("uninfected" = "#00BFC4", "infected" = "#F8766D")

# ---- Common theme (all bold, size 13) ----
base_theme <- theme_bw() +
  theme(
    strip.text   = element_text(face = "bold", size = 16, color = "black"),
    axis.text    = element_text(face = "bold", size = 16, color = "black"),
    axis.title   = element_text(face = "bold", size = 16, color = "black"),
    legend.title = element_text(face = "bold", size = 16, color = "black"),
    legend.text  = element_text(face = "bold", size = 16, color = "black"),
    legend.key.size = unit(0.5, "cm")
  )

# ---- Reusable panel builders ----------------------------------
# Top panel of a pair: x-axis removed (shared with bottom panel)
make_top <- function(data, yvar, ylab, tag) {
  ggplot(data, aes(x = Group, y = .data[[yvar]], fill = Status)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
                position = position_jitterdodge(jitter.width = 0.05,
                                                dodge.width = 0.8)) +
    scale_fill_manual(values = col_status) +
    scale_color_manual(values = col_status) +
    base_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.tag     = element_text(size = 16, face = "bold")) +
    labs(y = ylab, x = NULL,
         fill = "Infection status", color = "Infection status", tag = tag)
}

# Bottom panel of a pair: keeps the x-axis
make_bottom <- function(data, yvar, ylab, tag) {
  ggplot(data, aes(x = Group, y = .data[[yvar]], fill = Status)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
                position = position_jitterdodge(jitter.width = 0.05,
                                                dodge.width = 0.8)) +
    scale_fill_manual(values = col_status) +
    scale_color_manual(values = col_status) +
    base_theme +
    theme(plot.tag = element_text(size = 16, face = "bold")) +
    labs(y = ylab, x = NULL,
         fill = "Infection status", color = "Infection status", tag = tag)
}

# Single (standalone) panel: keeps the x-axis, no tag
make_single <- function(data, yvar, ylab) {
  ggplot(data, aes(x = Group, y = .data[[yvar]], fill = Status)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
                position = position_jitterdodge(jitter.width = 0.05,
                                                dodge.width = 0.8)) +
    scale_fill_manual(values = col_status) +
    scale_color_manual(values = col_status) +
    base_theme +
    labs(y = ylab, x = NULL,
         fill = "Infection status", color = "Infection status")
}

# ================================================================
# Figure 1: (a) mean aggregation proportion, (b) total aggregation time
# ================================================================
f1_a <- make_top(move_total, "mean_A",
                 "Mean aggregation\nproportion (%)", "a")
f1_b <- make_bottom(move_total, "total_time_agg",
                    "Total aggregation\ntime (min)", "b")

fig1 <- wrap_plots(f1_a, f1_b, ncol = 1) + plot_layout(guides = "collect")
fig1
ggsave("outputs/figures/Figure1_aggregation.pdf",
       plot = fig1, width = 8, height = 9)

# ================================================================
# Figure 2: (a) time to first aggregation, (b) mean total distance
#   Panel a uses only arenas where aggregation occurred
# ================================================================
f2_a <- make_top(filter(move_total, Aggregated == "yes"), "Start_Time",
                 "Time to first\naggregation (min)", "a")
f2_b <- make_bottom(move_total, "mean_gross_mm",
                    "Mean total\ndistance (mm)", "b")

fig2 <- wrap_plots(f2_a, f2_b, ncol = 1) + plot_layout(guides = "collect")
fig2
ggsave("outputs/figures/Figure2_time_distance.pdf",
       plot = fig2, width = 8, height = 9)

# ================================================================
# Figure S1: maximum aggregation proportion (single panel)
# ================================================================
figS1 <- make_single(move_total, "Amax",
                     "Maximum aggregation\nproportion (%)")
figS1
ggsave("outputs/figures/FigureS1_max_aggregation.pdf",
       plot = figS1, width = 8, height = 6)

# ================================================================
# Figure S2: (a) mean net displacement, (b) mean confinement index
# ================================================================
fS2_a <- make_top(move_total, "mean_net_mm",
                  "Mean net\ndisplacement (mm)", "a")
fS2_b <- make_bottom(move_total, "mean_CI",
                     "Mean confinement\nindex", "b")

figS2 <- wrap_plots(fS2_a, fS2_b, ncol = 1) + plot_layout(guides = "collect")
figS2
ggsave("outputs/figures/FigureS2_net_confinement.pdf",
       plot = figS2, width = 8, height = 9)

# ================================================================
# Figure S3: byssus thread count (single panel)
# ================================================================
figS3 <- make_single(move_total, "Byssus",
                     "Byssus thread count")
figS3
ggsave("outputs/figures/FigureS3_byssus.pdf",
       plot = figS3, width = 8, height = 6)

# ================================================================
# Figure S4: Overall Spearman correlations (n = 48 arenas)
#   (a) maximum aggregation proportion vs mean total distance
#   (b) maximum aggregation proportion vs mean net displacement
#   (c) maximum aggregation proportion vs byssus thread count
# ================================================================

# ---- Data: remove rows with missing values ----
vars_clean <- move_total %>%
  select(Amax, mean_gross_mm, mean_net_mm, Byssus) %>%
  drop_na()

# ---- Spearman correlations ----
cor_gross_overall  <- cor.test(vars_clean$Amax, vars_clean$mean_gross_mm,
                               method = "spearman", exact = FALSE)
cor_net_overall    <- cor.test(vars_clean$Amax, vars_clean$mean_net_mm,
                               method = "spearman", exact = FALSE)
cor_byssus_overall <- cor.test(vars_clean$Amax, vars_clean$Byssus,
                               method = "spearman", exact = FALSE)

# ---- Label helper: rho and p ----
make_label <- function(cor_obj) {
  paste0("rho = ", round(cor_obj$estimate, 2),
         ", p = ", signif(cor_obj$p.value, 2))
}

# ---- Correlation panel builder (same base_theme, size 16) ----
# left = TRUE keeps the y-axis title; FALSE drops it to avoid
# repeating the same label across the three horizontal panels
make_cor_panel <- function(data, xvar, xlab, cor_obj, tag, left = FALSE) {
  p <- ggplot(data, aes(x = .data[[xvar]], y = Amax)) +
    geom_point(size = 2.5, alpha = 0.5, color = "grey40") +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    annotate("text",
             x = min(data[[xvar]], na.rm = TRUE),
             y = max(data$Amax, na.rm = TRUE),
             label = make_label(cor_obj),
             hjust = 0, vjust = 1.2, size = 5) +
    base_theme +
    theme(plot.tag = element_text(size = 16, face = "bold")) +
    labs(x = xlab,
         y = "Maximum aggregation proportion (%)",
         tag = tag)
  if (!left) p <- p + theme(axis.title.y = element_blank())
  p
}

# ---- Three panels (only the first keeps the y-axis title) ----
s4_a <- make_cor_panel(vars_clean, "mean_gross_mm",
                       "Mean total distance (mm)",
                       cor_gross_overall, "a", left = TRUE)
s4_b <- make_cor_panel(vars_clean, "mean_net_mm",
                       "Mean net displacement (mm)",
                       cor_net_overall, "b", left = FALSE)
s4_c <- make_cor_panel(vars_clean, "Byssus",
                       "Byssus thread count",
                       cor_byssus_overall, "c", left = FALSE)

figS4 <- wrap_plots(s4_a, s4_b, s4_c, ncol = 3)
figS4
ggsave("outputs/figures/FigureS4_correlations.pdf",
       plot = figS4, width = 15, height = 6)
