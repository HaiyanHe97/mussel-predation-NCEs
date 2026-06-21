# Valve Gape Processing - 5-Minute Interval Averages
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   This script processes raw valve gape sensor readings into 5-min
#   interval averages for statistical analysis:
#     1) Raw gape values are converted to VC = sqrt(1 / Gape)
#     2) VC is normalised per individual mussel using 1-99% quantile
#        trimming followed by min-max scaling to 0-1 (VC_norm), to remove
#        sensor-specific extreme values and allow comparison across
#        individuals
#     3) The first 30 min of recording per mussel is split into six 5-min
#        bins, and VC_norm is averaged within each bin
#   Output: gap_5min_activity.csv, the dataset used as input for the valve
#   gape statistical models

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

View(total)

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

# 4. Calculate 5-minute gap activity per mussel
#    Step A: keep the first 30 min of recording per mussel (or actual
#    maximum if fewer data points are available)

total_first_30min <- total_mussel_id %>%
  group_by(Mussel_ID) %>%
  arrange(time) %>%
  slice(1:180) %>%   # if <180 rows, automatically keeps available rows
  mutate(
    time_start = min(time),
    minute_rel = as.numeric(difftime(time, time_start, units = "secs")) / 60,
    time_bin   = floor(minute_rel / 5) + 1   # bins 1-6
  ) %>%
  ungroup()

#    Step B: calculate mean VC_norm for each 5-min bin

gap_5min <- total_first_30min %>%
  group_by(Mussel_ID, time_bin) %>%
  summarise(
    VC_5min_avg = mean(VC_norm, na.rm = TRUE),
    n_points    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    time_mid_min = time_bin * 5
  )

#    Step C: retrieve metadata for each Mussel_ID

metadata <- total_mussel_id %>%
  group_by(Mussel_ID) %>%
  slice(1) %>%
  ungroup()

#    Step D: merge

gap_5min_full <- gap_5min %>%
  left_join(metadata, by = "Mussel_ID") %>%
  select(Mussel_ID, Group, Status, Rep, time_bin, time_mid_min, VC_5min_avg,
         n_points, everything())

View(gap_5min_full)

# 5. Export results

write.csv(
  gap_5min_full,
  "data/processed/gap_5min_activity.csv",
  row.names = FALSE
)
