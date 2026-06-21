# Mussel Movement & Aggregation - Figures
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Plot order follows the analysis code order
# p1-p8: main boxplots (Group x Status)
# p3b, p4b: alternative plots for total_time_agg and mean_A (Status as x-axis)
# p3c, p4c: alternative plots for total_time_agg and mean_A (Group as x-axis,
#           Status as jitter colour)
# p_cor: overall Spearman correlation figure (main)
# p_cor_supp: grouped correlation figure (supplementary)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(purrr)

# 0. Data preparation

move_total <- read.csv(
  "data/processed/merged_movement.csv",
  header = TRUE, sep = ",", stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

move_total$Group  <- factor(move_total$Group,
                            levels = c("Control", "CM", "HT", "HS"))
move_total$Status <- factor(move_total$Status,
                            levels = c("uninfected", "infected"))

# Colour palettes 
col_status <- c("uninfected" = "#00BFC4", "infected" = "#F8766D")
col_group  <- c("Control" = "#1f78b4", "CM" = "#33a02c",
                "HT" = "#e31a1c", "HS" = "#ff7f00")

# Common theme 
base_theme <- theme_bw() +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14, color = "black"),
    axis.title   = element_text(face = "bold", size = 18, color = "black"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text  = element_text(face = "bold", size = 16, color = "black")
  )

# p1: Time to first aggregation - Start_Time (min)
#     Only arenas where aggregation was observed (Aggregated == "yes")

p1 <- move_total %>%
  filter(Aggregated == "yes") %>%
  ggplot(aes(x = Group, y = Start_Time, fill = Status)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 1,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    alpha = 0.5, size = 1.5,
    position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = col_status) +
  base_theme +
  labs(y = "Time to first aggregation (min)", x = " ",
       fill = "Infection status")

p1

ggsave(
  "outputs/figures/p1_start_time.pdf",
  plot = p1, width = 8, height = 6
)

# p2: Mean total distance - mean_gross_mm (mm)

p2 <- ggplot(move_total,
             aes(x = Group, y = mean_gross_mm, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  base_theme +
  labs(y = "Mean total distance (mm)", x = " ",
       fill = "Infection status")

p2

ggsave(
  "outputs/figures/p2_total_distance.pdf",
  plot = p2, width = 8, height = 6
)

# p3: Total aggregation time - total_time_agg (min)
#     Version A: Group as x-axis, Status as fill

p3 <- ggplot(move_total,
             aes(x = Group, y = total_time_agg, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  scale_color_manual(values = col_status) +
  base_theme +
  labs(y = "Total aggregation time (min)", x = " ",
       fill = "Infection status", color = "Infection status")

p3

ggsave(
  "outputs/figures/p3_total_agg_time.pdf",
  plot = p3, width = 8, height = 6
)

# p3b: Total aggregation time - total_time_agg (min)
#      Version B: Status as x-axis, Group as colour
#      Rationale: highlights the infection status effect independently
#                 from cue type

p3b <- ggplot(move_total,
              aes(x = Status, y = total_time_agg)) +
  geom_boxplot(aes(group = Status), fill = "grey80", color = "black",
               outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.15, alpha = 0.5, size = 2) +
  scale_color_manual(values = col_group) +
  base_theme +
  labs(y = "Total aggregation time (min)", x = " ",
       color = "Predation risk cue")

p3b

ggsave(
  "outputs/figures/p3b_total_agg_time_status.pdf",
  plot = p3b, width = 6, height = 6
)

# p3c: Total aggregation time - total_time_agg (min)
#      Group as x-axis, grey boxplot, Status as jitter colour
#      Rationale: highlights cue type effect with infection status shown as points

p3c <- ggplot(move_total,
              aes(x = Group, y = total_time_agg)) +
  geom_boxplot(aes(group = Group), fill = "grey80", color = "black",
               outlier.shape = NA) +
  geom_jitter(aes(color = Status), width = 0.15, alpha = 0.5, size = 2) +
  scale_color_manual(values = col_status) +
  base_theme +
  labs(y = "Total aggregation time (min)", x = " ",
       color = "Infection status")

p3c

ggsave(
  "outputs/figures/p3c_total_agg_time_group.pdf",
  plot = p3c, width = 6, height = 6
)

# p4: Mean aggregation proportion - mean_A (%)
#     Version A: Group as x-axis, Status as fill

p4 <- ggplot(move_total,
             aes(x = Group, y = mean_A, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  scale_color_manual(values = col_status) +
  base_theme +
  labs(y = "Mean aggregation proportion (%)", x = " ",
       fill = "Infection status", color = "Infection status")

p4

ggsave(
  "outputs/figures/p4_mean_A.pdf",
  plot = p4, width = 8, height = 6
)

# p4b: Mean aggregation proportion - mean_A (%)
#      Version B: Status as x-axis, Group as colour
#      Rationale: highlights the infection status effect independently
#                 from cue type

p4b <- ggplot(move_total,
              aes(x = Status, y = mean_A)) +
  geom_boxplot(aes(group = Status), fill = "grey80", color = "black",
               outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.15, alpha = 0.5, size = 2) +
  scale_color_manual(values = col_group) +
  base_theme +
  labs(y = "Mean aggregation proportion (%)", x = " ",
       color = "Predation risk cue")

p4b

ggsave(
  "outputs/figures/p4b_mean_A_status.pdf",
  plot = p4b, width = 6, height = 6
)

# p4c: Mean aggregation proportion - mean_A (%)
#      Group as x-axis, grey boxplot, Status as jitter colour
#      Rationale: highlights cue type effect with infection status shown as points

p4c <- ggplot(move_total,
              aes(x = Group, y = mean_A)) +
  geom_boxplot(aes(group = Group), fill = "grey80", color = "black",
               outlier.shape = NA) +
  geom_jitter(aes(color = Status), width = 0.15, alpha = 0.5, size = 2) +
  scale_color_manual(values = col_status) +
  base_theme +
  labs(y = "Mean aggregation proportion (%)", x = " ",
       color = "Infection status")

p4c

ggsave(
  "outputs/figures/p4c_mean_A_group.pdf",
  plot = p4c, width = 6, height = 6
)

# p5: Byssus thread count - Byssus (count per mussel)

p5 <- ggplot(move_total,
             aes(x = Group, y = Byssus, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  base_theme +
  labs(y = "Byssus thread", x = " ",
       fill = "Infection status")

p5

ggsave(
  "outputs/figures/p5_byssus.pdf",
  plot = p5, width = 8, height = 6
)

# p6: Mean net displacement - mean_net_mm (mm)

p6 <- ggplot(move_total,
             aes(x = Group, y = mean_net_mm, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  base_theme +
  labs(y = "Mean net displacement (mm)", x = " ",
       fill = "Infection status")

p6

ggsave(
  "outputs/figures/p6_net_distance.pdf",
  plot = p6, width = 8, height = 6
)

# p7: Mean confinement index - mean_CI (ratio, 0-1)

p7 <- ggplot(move_total,
             aes(x = Group, y = mean_CI, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  base_theme +
  labs(y = "Mean confinement index", x = " ",
       fill = "Infection status")

p7

ggsave(
  "outputs/figures/p7_confinement_index.pdf",
  plot = p7, width = 8, height = 6
)

# p8: Maximum aggregation proportion - Amax (%)

p8 <- ggplot(move_total,
             aes(x = Group, y = Amax, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 1,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Status), alpha = 0.5, size = 1.5,
              position = position_jitterdodge(jitter.width = 0.05,
                                              dodge.width = 0.8)) +
  scale_fill_manual(values = col_status) +
  scale_color_manual(values = col_status) +
  base_theme +
  labs(y = "Maximum aggregation proportion (%)", x = " ",
       fill = "Infection status", color = "Infection status")

p8

ggsave(
  "outputs/figures/p8_Amax.pdf",
  plot = p8, width = 8, height = 6
)

# p_cor: Overall Spearman correlation figure (MAIN)
#        n = 48 arenas across all treatments
#        Panel A: Amax ~ mean total distance (taxis)
#        Panel B: Amax ~ mean net displacement (taxis)
#        Panel C: Amax ~ byssus thread count (aggregation strength)

vars_clean <- move_total %>%
  select(Amax, mean_gross_mm, mean_net_mm, Byssus) %>%
  drop_na()

cor_gross_overall  <- cor.test(vars_clean$Amax, vars_clean$mean_gross_mm,
                               method = "spearman", exact = FALSE)
cor_net_overall    <- cor.test(vars_clean$Amax, vars_clean$mean_net_mm,
                               method = "spearman", exact = FALSE)
cor_byssus_overall <- cor.test(vars_clean$Amax, vars_clean$Byssus,
                               method = "spearman", exact = FALSE)

make_label <- function(cor_obj) {
  paste0("\u03c1 = ", round(cor_obj$estimate, 2),
         ", p = ", signif(cor_obj$p.value, 2))
}

cor_theme <- theme_bw() +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14, color = "black"),
    axis.title   = element_text(face = "bold", size = 18, color = "black"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text  = element_text(face = "bold", size = 16, color = "black"),
    legend.position = "bottom"
  )

# Panel A: Amax vs mean total distance
pa <- ggplot(vars_clean, aes(x = mean_gross_mm, y = Amax)) +
  geom_point(size = 2.5, alpha = 0.5, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = min(vars_clean$mean_gross_mm, na.rm = TRUE),
           y = max(vars_clean$Amax, na.rm = TRUE),
           label = make_label(cor_gross_overall),
           hjust = 0, vjust = 1.2, size = 5) +
  cor_theme +
  labs(x = "Mean total distance (mm)",
       y = "Maximum aggregation proportion (%)")

# Panel B: Amax vs mean net displacement
pb <- ggplot(vars_clean, aes(x = mean_net_mm, y = Amax)) +
  geom_point(size = 2.5, alpha = 0.5, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = min(vars_clean$mean_net_mm, na.rm = TRUE),
           y = max(vars_clean$Amax, na.rm = TRUE),
           label = make_label(cor_net_overall),
           hjust = 0, vjust = 1.2, size = 5) +
  cor_theme +
  labs(x = "Mean net displacement (mm)",
       y = "Maximum aggregation proportion (%)")

# Panel C: Amax vs byssus thread count
pc <- ggplot(vars_clean, aes(x = Byssus, y = Amax)) +
  geom_point(size = 2.5, alpha = 0.5, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = min(vars_clean$Byssus, na.rm = TRUE),
           y = max(vars_clean$Amax, na.rm = TRUE),
           label = make_label(cor_byssus_overall),
           hjust = 0, vjust = 1.2, size = 5) +
  cor_theme +
  labs(x = "Byssus thread",
       y = "Maximum aggregation proportion (%)")

p_cor <- ggarrange(pa, pb, pc,
                   ncol = 3, nrow = 1,
                   labels = c("A", "B", "C"))
p_cor

ggsave(
  "outputs/figures/p_cor_overall.pdf",
  plot = p_cor, width = 15, height = 8
)
ggsave(
  "outputs/figures/p_cor_overall.png",
  plot = p_cor, width = 15, height = 7,
  dpi = 300
)

# p_cor_supp: Grouped Spearman correlation figure (SUPPLEMENTARY)
#             n = 6 per Group x Status subgroup - exploratory only
#             Faceted by Status (rows) x Variable (columns)
#             Coloured by Group

plot_data_supp <- move_total %>%
  select(Group, Status, Amax, mean_net_mm, mean_gross_mm, Byssus) %>%
  pivot_longer(
    cols = c(mean_gross_mm, mean_net_mm, Byssus),
    names_to  = "Variable",
    values_to = "Value"
  ) %>%
  mutate(Variable = factor(Variable,
                           levels = c("mean_gross_mm", "mean_net_mm", "Byssus"),
                           labels = c("Mean total distance (mm)",
                                      "Mean net displacement (mm)",
                                      "Byssus thread count")))

p_cor_supp <- ggplot(plot_data_supp,
                     aes(x = Value, y = Amax, color = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(Status ~ Variable, scales = "free_x") +
  scale_color_manual(values = col_group) +
  theme_bw() +
  theme(
    legend.position  = "top",
    strip.text       = element_text(face = "bold", size = 16),
    axis.text        = element_text(face = "bold", size = 14, color = "black"),
    axis.title       = element_text(face = "bold", size = 18, color = "black"),
    legend.title     = element_text(face = "bold", size = 14, color = "black"),
    legend.text      = element_text(face = "bold", size = 16, color = "black")
  ) +
  labs(
    x       = " ",
    y       = "Maximum aggregation proportion (%)",
    color   = "Predation risk cue",
    caption = "Exploratory analysis: n = 6 per subgroup. Interpret with caution."
  )

p_cor_supp

ggsave(
  "outputs/figures/p_cor_supp.pdf",
  plot = p_cor_supp, width = 12, height = 7
)
