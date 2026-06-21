# Heart Rate Data Processing - Merge with Metadata
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description: Merge raw heartbeatr output (one file per batch, generated
#              by 01_calc_heart_rate.R) with experimental metadata, clean
#              duplicates and known recording errors, convert Hz to bpm,
#              and prepare the final dataset for analysis.

library(dplyr)
library(readr)

# Define paths

path    <- "data/processed/heart_rate_by_batch"
raw_csv <- "data/raw/heart_rate_06-07_raw.csv"

# 1. Read raw heartbeatr output
# Only read files ending in _heartrate.csv to exclude previously merged files

csv_files <- list.files(path, pattern = "_heartrate\\.csv$", full.names = TRUE)
cat("Number of heartrate files found:", length(csv_files), "\n")  # should be 29

merged_df <- csv_files %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows()

# Remove index column if present
if ("...1" %in% names(merged_df)) {
  merged_df <- merged_df %>% select(-"...1")
}

cat("Rows in merged_df:", nrow(merged_df), "\n")

# 2. Standardise id format
# Capitalise first letter to match df_raw (e.g. "mussel_1" -> "Mussel_1")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

merged_df <- merged_df %>%
  mutate(id = capitalize_first(id))

# 2b. Correct known id recording error
# Mussel_286 (Rep 29, HT uninfected, 2025-07-01 09:15) was incorrectly
# recorded as Mussel_256 in the heartrate output file. Corrected here
# based on experiment date and time.

merged_df <- merged_df %>%
  mutate(id = case_when(
    id == "Mussel_256" & as.Date(time) == as.Date("2025-07-01") ~ "Mussel_286",
    TRUE ~ id
  ))

# Verify correction
cat("Mussel_256 records remaining:", nrow(subset(merged_df, id == "Mussel_256")), "\n")  # should be 5-8
cat("Mussel_286 records created:",   nrow(subset(merged_df, id == "Mussel_286")), "\n")  # should be 5-8

# 3. Read experimental metadata

df_raw <- read_delim(raw_csv, delim = ";", show_col_types = FALSE) %>%
  mutate(id = trimws(id))

cat("Rows in df_raw:", nrow(df_raw), "\n")

# 4. Join metadata to heartrate data
# Left join: keep all heartrate rows, attach metadata where id matches

merged_final <- merged_df %>%
  left_join(df_raw, by = "id")

# Check how many rows have no metadata match
n_na <- sum(is.na(merged_final$Rep))
cat("Rows with no Rep match (will be removed):", n_na, "\n")

# Remove rows with no metadata match
merged_final <- merged_final %>%
  filter(!is.na(Rep))

# 5. Convert Hz to bpm

merged_final <- merged_final %>%
  mutate(bpm = hz * 60)

# 6. Check for duplicates

n_dupes <- merged_final %>%
  group_by(id, i) %>%
  filter(n() > 1) %>%
  nrow()
cat("Duplicate id x i rows remaining:", n_dupes, "\n")  # should be 0

# 7. Check sample sizes

cat("\nSample sizes (individuals per Group x Status):\n")
with(unique(merged_final[, c("id", "Group", "Status")]), table(Group, Status))

# 8. Save cleaned file

write_csv(merged_final, "data/processed/heartrate_merged_final.csv")
cat("\nSaved: heartrate_merged_final.csv\n")
