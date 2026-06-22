# Mussel Movement & Aggregation Analysis
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Design: 2 (Status: infected / uninfected) x 4 (Group: Control / CM / HT / HS)
# Total: 8 treatments x 6 replicates = 48 arenas, 7 mussels per arena
# Observation unit: arena (petri dish)
# Random effect: Rep (experimental block, 6 levels)

library(dplyr)
library(lme4)        # linear mixed models (LMM)
library(lmerTest)    # p-values for LMM via Satterthwaite approximation
library(emmeans)     # estimated marginal means & post-hoc comparisons
library(car)         # Levene's test
library(ARTool)      # Aligned Rank Transform non-parametric two-way ANOVA
library(DHARMa)      # simulation-based residual diagnostics for mixed models
library(tidyr)       # nest/unnest for grouped correlations
library(purrr)       # map() for grouped correlations

# Read & format data

move_total <- read.csv(
  "data/processed/merged_movement.csv",
  header = TRUE, sep = ",", stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

# Factor encoding
move_total$Group  <- factor(move_total$Group,
                            levels = c("Control", "CM", "HT", "HS"))
move_total$Status <- factor(move_total$Status,
                            levels = c("uninfected", "infected"))
move_total$Rep    <- factor(move_total$Rep)

# Numeric coercion
move_total$sum_intensity  <- as.numeric(move_total$sum_intensity)
move_total$n_infected     <- as.numeric(move_total$n_infected)
move_total$Byssus         <- as.numeric(move_total$Byssus)
move_total$Start_Time     <- as.numeric(move_total$Start_Time)
move_total$Amax           <- as.numeric(move_total$Amax)
move_total$mean_A         <- as.numeric(move_total$mean_A)
move_total$total_time_agg <- as.numeric(move_total$total_time_agg)

# 1. Time to First Aggregation - Start_Time (min)
#    Definition: latency from trial start until the first 5-min interval at
#                which any aggregation was observed (proportion > 0).
#                Arenas where no aggregation occurred throughout the trial
#                have NA and are excluded from this analysis.
#    Model: LMM with log(x+1) transformation
#    Rationale: right-skewed distribution with structural zeros excluded;
#               log-transform improves residual normality.
#               ART not used here because structural NAs would distort
#               rank alignment.
#    Random effect: Rep (experimental block)

move_total_clean <- move_total[!is.na(move_total$Start_Time), ]

lmm_time <- lmer(log(Start_Time + 1) ~ Status * Group + (1 | Rep),
                 data = move_total_clean)
summary(lmm_time)
anova(lmm_time)
plot(lmm_time)

sim_res <- simulateResiduals(fittedModel = lmm_time)
plot(sim_res)

# Post-hoc: infected vs uninfected within each cue type
emm9  <- emmeans(lmm_time, ~ Status | Group)
pairs(emm9)
# Post-hoc: cue type effect within each infection status
emm10 <- emmeans(lmm_time, ~ Group | Status)
pairs(emm10)


# 2. Total Aggregation Time - total_time_agg (min)
#    Definition: cumulative duration (min) during which aggregation was
#                observed, calculated as the number of 5-min intervals with
#                a non-zero aggregation proportion multiplied by 5.
#                A non-zero proportion indicates at least two mussels were
#                in direct physical contact during that interval.
#    Model: ART (Aligned Rank Transform) two-way ANOVA
#    Rationale: discrete values (multiples of 5), zero-inflated, non-normal
#               distribution confirmed by Shapiro-Wilk and Levene tests below.
#               ART is the non-parametric equivalent of two-way ANOVA and
#               preserves the factorial design (Group x Status).
#    Note: interaction term is exploratory - only interpret if significant

model_art <- art(total_time_agg ~ Group * Status, data = move_total)
anova(model_art)

# Normality and homogeneity checks to justify use of ART
shapiro.test(residuals(model_art))
leveneTest(total_time_agg ~ Group * Status, data = move_total)

# Post-hoc: cue type effect
model_group <- artlm(model_art, "Group")
emmeans(model_group, pairwise ~ Group)

# Post-hoc: infection status effect
model_status <- artlm(model_art, "Status")
emmeans(model_status, pairwise ~ Status)

# Post-hoc: interaction (uncomment only if interaction term is significant)
# model_inter <- artlm(model_art, "Group:Status")
# emmeans(model_inter, pairwise ~ Group | Status)
# emmeans(model_inter, pairwise ~ Status | Group)

# 3. Mean Aggregation Proportion - mean_A (%)
#    Definition: arithmetic mean of per-interval aggregation proportions
#                across all 5-min observation intervals in the trial.
#                Reflects the average tendency to aggregate throughout
#                the full duration of the experiment.
#    Model: ART two-way ANOVA
#    Rationale: proportion data (0-100%), non-normal distribution.
#               ART preferred over LMM on raw or transformed proportions.
#    Note: interaction term fitted but only interpret if significant

model_art <- art(mean_A ~ Group * Status, data = move_total)
anova(model_art)

# Post-hoc: cue type effect
model_group <- artlm(model_art, "Group")
emmeans(model_group, pairwise ~ Group)

# Post-hoc: infection status effect
model_status <- artlm(model_art, "Status")
emmeans(model_status, pairwise ~ Status)

# Post-hoc: interaction (uncomment only if interaction term is significant)
# model_inter <- artlm(model_art, "Group:Status")
# emmeans(model_inter, pairwise ~ Group | Status)
# emmeans(model_inter, pairwise ~ Status | Group)

# 4. Maximum Aggregation Proportion
#    Definition: the highest aggregation proportion recorded across all 5-min
#                observation intervals during the trial, representing the peak
#                collective aggregation response. Also used as the response
#                variable in conspecific taxis and aggregation strength
#                correlations (Section 10).
#    Model: ART two-way ANOVA
#    Rationale: proportion data (0-100%), non-normal distribution.

model_art <- art(Amax ~ Group * Status, data = move_total)
anova(model_art)

# Post-hoc: cue type effect
model_group <- artlm(model_art, "Group")
emmeans(model_group, pairwise ~ Group)

# Post-hoc: infection status effect
model_status <- artlm(model_art, "Status")
emmeans(model_status, pairwise ~ Status)

# Post-hoc: interaction (uncomment only if interaction term is significant)
# model_inter <- artlm(model_art, "Group:Status")
# emmeans(model_inter, pairwise ~ Group | Status)

# 5. Mean Total Distance - mean_gross_mm (mm)
#    Definition: mean total path length (mm) travelled by each mussel across
#                the full trial, summing all positional steps regardless of
#                direction. Averaged across the 7 individuals per arena.
#                Reflects overall locomotor activity.
#    Model: LMM (untransformed - residuals checked below)
#    Random effect: Rep (experimental block)

lmm_gross <- lmer(mean_gross_mm ~ Group * Status + (1 | Rep),
                  data = move_total)
summary(lmm_gross)
anova(lmm_gross)

sim_res <- simulateResiduals(fittedModel = lmm_gross)
plot(sim_res)

# Post-hoc: infected vs uninfected within each cue type
emm3 <- emmeans(lmm_gross, ~ Status | Group)
pairs(emm3)
# Post-hoc: cue type effect within each infection status
emm4 <- emmeans(lmm_gross, ~ Group | Status)
pairs(emm4)


# 6. Mean Net Displacement - mean_net_mm (mm)
#    Definition: mean straight-line distance (mm) between each mussel's
#                starting and ending position, averaged across the 7 individuals
#                per arena. Unlike total distance, this metric captures net
#                directional relocation: a mussel that moves far but returns
#                to its start has high total distance but low net displacement.
#    Model: LMM with log(x+1) transformation
#    Rationale: right-skewed with zeros; log-transform improves normality.
#    Random effect: Rep (experimental block)

lmm_net2 <- lmer(log(mean_net_mm + 1) ~ Group * Status + (1 | Rep),
                 data = move_total)
summary(lmm_net2)
anova(lmm_net2)
plot(lmm_net2)

sim_res <- simulateResiduals(fittedModel = lmm_net2)
plot(sim_res)

# Post-hoc: infected vs uninfected within each cue type
emm5 <- emmeans(lmm_net2, ~ Status | Group)
pairs(emm5)
# Post-hoc: cue type effect within each infection status
emm6 <- emmeans(lmm_net2, ~ Group | Status)
pairs(emm6)

# 7. Mean Confinement Index - mean_CI (ratio, 0-1)
#    Definition: ratio of net displacement to total distance for each mussel,
#                averaged across the 7 individuals per arena.
#                Values near 1 = straight directed movement.
#                Values near 0 = tortuous or local movement with little net
#                progress (e.g. circling or remaining in place).
#                Higher values suggest stronger directional movement toward
#                conspecifics or in response to predator cues.
#    Model: LMM (untransformed - residuals checked below)
#    Random effect: Rep (experimental block)

lmm_ci <- lmer(mean_CI ~ Group * Status + (1 | Rep),
               data = move_total)
summary(lmm_ci)
anova(lmm_ci)
plot(lmm_ci)
qqnorm(resid(lmm_ci)); qqline(resid(lmm_ci))

sim_res <- simulateResiduals(fittedModel = lmm_ci)
plot(sim_res)

# Post-hoc: infected vs uninfected within each cue type
emm7 <- emmeans(lmm_ci, ~ Status | Group)
pairs(emm7)
# Post-hoc: cue type effect within each infection status
emm8 <- emmeans(lmm_ci, ~ Group | Status)
pairs(emm8)

# 8. Byssus Thread Count - Byssus (count per mussel)
#    Definition: number of byssus threads attached to the substratum or to
#                other mussels per individual, counted at the end of the trial.
#                Used as a proxy for settlement tendency and commitment to
#                remaining in place (aggregation reinforcement).
#    Model: LMM with log(x+1) transformation
#    Rationale: count data, right-skewed with zeros; log-transform improves
#               residual normality.
#    Random effect: Rep (experimental block)

lmm_by2 <- lmer(log(Byssus + 1) ~ Group * Status + (1 | Rep),
                data = move_total)
summary(lmm_by2)
anova(lmm_by2)

sim_res <- simulateResiduals(fittedModel = lmm_by2)
plot(sim_res)

# Post-hoc: infected vs uninfected within each cue type (Bonferroni: conservative)
emm1 <- emmeans(lmm_by2, ~ Status | Group, adjust = "bonferroni")
pairs(emm1)
# Post-hoc: cue type effect within each infection status (Tukey)
emm2 <- emmeans(lmm_by2, ~ Group | Status, adjust = "tukey")
pairs(emm2)







# 9. Conspecific Taxis & Aggregation Strength - Spearman correlations
#
#     (A) Conspecific taxis: does higher locomotion predict higher aggregation?
#         Variables: Amax ~ mean_gross_mm / mean_net_mm
#         Expected: positive correlation - mussels that travel further are more
#         likely to encounter and join conspecifics, resulting in higher Amax.
#
#     (B) Aggregation strength: does higher aggregation predict more byssus?
#         Variables: Amax ~ Byssus
#         Expected: positive correlation - more aggregated mussels produce more
#         byssus threads, reflecting stronger collective binding force.
#
#     Method: Spearman rank correlation (non-parametric; consistent with ART
#             used for main analyses; exact = FALSE to handle ties in
#             proportion and count data).
#
#     Two levels of analysis:
#     - Overall (n = 48 arenas): primary analysis
#     - Within each Group x Status combination (n = 6 per subgroup): exploratory
#       only, due to limited statistical power at small sample sizes

# Part A: Overall Spearman correlations (n = 48)

# Remove rows with missing values in relevant variables
vars_clean <- move_total %>%
  select(Amax, mean_gross_mm, mean_net_mm, Byssus) %>%
  drop_na()

# Taxis: Amax vs mean total distance
cor_gross_overall <- cor.test(vars_clean$Amax,
                              vars_clean$mean_gross_mm,
                              method = "spearman",
                              exact = FALSE)
cor_gross_overall

# Taxis: Amax vs mean net displacement
cor_net_overall <- cor.test(vars_clean$Amax,
                            vars_clean$mean_net_mm,
                            method = "spearman",
                            exact = FALSE)
cor_net_overall

# Aggregation strength: Amax vs byssus thread count
cor_byssus_overall <- cor.test(vars_clean$Amax,
                               vars_clean$Byssus,
                               method = "spearman",
                               exact = FALSE)
cor_byssus_overall

# Part B: Spearman correlations within each Group x Status (exploratory)

# Helper function: run Spearman correlation and return rho and p-value
run_cor <- function(data, x, y) {
  test <- cor.test(data[[x]], data[[y]], method = "spearman", exact = FALSE)
  data.frame(
    rho = as.numeric(test$estimate),
    p   = test$p.value
  )
}

# Taxis: Amax vs mean total distance
cor_gross <- move_total %>%
  group_by(Group, Status) %>%
  nest() %>%
  mutate(res = map(data, ~ run_cor(.x, "Amax", "mean_gross_mm"))) %>%
  unnest(res) %>%
  select(Group, Status, rho, p)

# Taxis: Amax vs mean net displacement
cor_net <- move_total %>%
  group_by(Group, Status) %>%
  nest() %>%
  mutate(res = map(data, ~ run_cor(.x, "Amax", "mean_net_mm"))) %>%
  unnest(res) %>%
  select(Group, Status, rho, p)

# Aggregation strength: Amax vs byssus thread count
cor_byssus <- move_total %>%
  group_by(Group, Status) %>%
  nest() %>%
  mutate(res = map(data, ~ run_cor(.x, "Amax", "Byssus"))) %>%
  unnest(res) %>%
  select(Group, Status, rho, p)

# Part C: Summary table with significance stars

add_sig <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
}

final_table <- cor_gross %>%
  rename(rho_gross = rho, p_gross = p) %>%
  left_join(cor_net    %>% rename(rho_net    = rho, p_net    = p),
            by = c("Group", "Status")) %>%
  left_join(cor_byssus %>% rename(rho_byssus = rho, p_byssus = p),
            by = c("Group", "Status")) %>%
  mutate(
    sig_gross  = add_sig(p_gross),
    sig_net    = add_sig(p_net),
    sig_byssus = add_sig(p_byssus)
  )

final_table

