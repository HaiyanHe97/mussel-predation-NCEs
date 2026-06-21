# Valve Gape Processing - 5-Minute Interval CV
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   This script processes raw valve gape sensor readings into 5-min
#   interval coefficients of variation (CV) for statistical analysis:
#     1) Raw gape values are converted to VC = sqrt(1 / Gape)
#     2) VC is normalised per individual mussel using 1-99% quantile
#        trimming followed by min-max scaling to 0-1 (VC_norm), to remove
#        sensor-specific extreme values and allow comparison across
#        individuals
#     3) The first 30 min of recording per mussel is split into six 5-min
#        bins, and the coefficient of variation (SD/mean) of VC_norm is
#        calculated within each bin to capture temporal variability in
#        valve opening behaviour
#   Output: gap_5min_CV.csv, the dataset used as input for the valve gape
#   CV statistical models

library(dplyr)
library(lubridate)
library(readr)

# 1. Load data

total <- read.csv(
  "data/processed/all_gape_with_mussel.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

# 2. Build timestamp and calculate VC

total <- total %>%
  mutate(
    time = as.POSIXct(
      paste(Year, Month, Day, Hour, Minute, Second, sep = "-"),
      format = "%Y-%m-%d-%H-%M-%S",
      tz = "UTC"
    )
  ) %>%
  mutate(VC = sqrt(1 / Gape))

total$Rep       <- factor(total$Rep)
total$Mussel_ID <- factor(total$Mussel_ID)
total$Group     <- factor(total$Group, levels = c("Control", "CM", "HT", "HS"))
total$Status    <- factor(total$Status, levels = c("uninfected", "infected"))

View(total)

# 3. VC normalisation by 1-99% quantile trimming
# Trims extreme values per individual to the 1st-99th percentile range,
# then rescales to 0-1 (VC_norm), to account for sensor-specific baseline
# differences between individuals

total_mussel_id <- total %>%
  group_by(Mussel_ID) %>%
  mutate(
    q_low   = quantile(VC, 0.01, na.rm = TRUE),
    q_high  = quantile(VC, 0.99, na.rm = TRUE),
    VC_trim = pmin(pmax(VC, q_low), q_high),
    VC_norm = (VC_trim - q_low) / (q_high - q_low)
  ) %>%
  ungroup()

# 4. Calculate 5-minute CV per mussel
#    Step A: keep the first 30 min of recording per mussel, split into
#    six 5-min bins (bins beyond 6 are discarded)

total_first_30min <- total_mussel_id %>%
  group_by(Mussel_ID) %>%
  arrange(time) %>%
  slice(1:180) %>%   # if <180 rows, automatically keeps available rows
  mutate(
    time_start = min(time),
    minute_rel = as.numeric(difftime(time, time_start, units = "secs")) / 60,
    time_bin   = floor(minute_rel / 5) + 1   # convert 0-5 to bins 1-6
  ) %>%
  ungroup() %>%
  filter(time_bin >= 1 & time_bin <= 6)   # remove any bin > 6

#    Step B: calculate CV (SD/mean) of VC_norm for each 5-min bin

gap_5min <- total_first_30min %>%
  group_by(Mussel_ID, time_bin) %>%
  summarise(
    n_points = n(),
    VC_5min_CV = ifelse(mean(VC_norm, na.rm = TRUE) != 0,
                        sd(VC_norm, na.rm = TRUE) / mean(VC_norm, na.rm = TRUE),
                        NA),
    .groups = "drop"
  ) %>%
  mutate(
    time_mid_min = time_bin * 5   # midpoint of interval for bins 1-6
  )

#    Step C: retrieve metadata for each Mussel_ID

metadata <- total_mussel_id %>%
  group_by(Mussel_ID) %>%
  slice(1) %>%
  ungroup()

#    Step D: merge

gap_5min_CV <- gap_5min %>%
  left_join(metadata, by = "Mussel_ID") %>%
  select(Mussel_ID, Group, Status, Rep, time_bin, time_mid_min, n_points, everything())

View(gap_5min_CV)

# 5. Export results

write.csv(
  gap_5min_CV,
  "data/processed/gap_5min_CV.csv",
  row.names = FALSE
)
