# =============================================================================
# Infection Intensity Summary
# Project: Effects of parasitic infection on predation risk responses
#          in blue mussels (Mytilus edulis)
#
# Purpose:
#   Summarise infection intensity (number of encysted metacercariae per mussel)
#   across all infected individuals used in the study, pooled over the three
#   behavioural assays (movement, valve gape, heart rate).
#
# Notes:
#   - Status    : infection status confirmed by dissection (filter variable)
#   - Intensity : number of metacercariae counted per mussel
#   - Mussel_ID : individual label (e.g. G25); shared across all three assays
#   - Mussels were reused across assays, so an individual may appear in more
#     than one dataset. We de-duplicate by Mussel_ID, keeping the first
#     non-missing Intensity, so each infected mussel is counted once.
#   - Heart-rate signal-quality exclusions are NOT applied here: they concern
#     cardiac data only and are unrelated to infection intensity.
# =============================================================================

# 1. Libraries ----------------------------------------------------------------
library(dplyr)

# 2. Load data ----------------------------------------------------------------
# Movement raw data is semicolon-separated; the processed files use commas.
move_raw <- read.csv("data/raw/movement_row_clean.csv",
                     sep = ";", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
View(move_raw)
valve_total <- read.csv("data/processed/gap_5min_activity.csv",
                        sep = ",", stringsAsFactors = FALSE, fileEncoding = "UTF-8")

hr_total <- read.csv("data/processed/heartrate_merged_final.csv",
                     sep = ",", stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 3. Pool infected individuals and de-duplicate by Mussel_ID ------------------
all_infected <- bind_rows(
    move_raw    %>% select(Mussel_ID, Status, Intensity),
    valve_total %>% select(Mussel_ID, Status, Intensity),
    hr_total    %>% select(Mussel_ID, Status, Intensity)
  ) %>%
  mutate(Intensity = as.numeric(Intensity)) %>%
  filter(Status == "infected") %>%
  group_by(Mussel_ID) %>%
  summarise(Intensity = first(na.omit(Intensity)), .groups = "drop")

# 4. Data check: individuals with no recorded intensity -----------------------
all_infected %>% summarise(n_total = n(), n_missing = sum(is.na(Intensity)))

# 5. Summary statistics -------------------------------------------------------
# Intensity is right-skewed; median + range are the reported values.
all_infected %>%
  filter(!is.na(Intensity)) %>%
  summarise(
    n      = n(),
    min    = min(Intensity),
    max    = max(Intensity),
    mean   = mean(Intensity),
    sd     = sd(Intensity),
    median = median(Intensity)
  )



all_infected %>% filter(!is.na(Intensity)) %>%
  summarise(median = median(Intensity),
            Q1 = quantile(Intensity, 0.25),
            Q3 = quantile(Intensity, 0.75))
