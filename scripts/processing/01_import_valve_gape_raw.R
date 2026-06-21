# Valve Gape Data Import - Raw Sensor Data to Long Format
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   This script reads raw multi-channel valve gape sensor recordings (one
#   CSV per sensor unit, A/B, per replicate), reshapes them from wide to
#   long format (one row per channel per timestamp), and matches each
#   channel to its corresponding mussel using metadata (Sensor_ID, Rep,
#   Line_ID). Rows with no matching Mussel_ID (i.e. unused/empty channels)
#   are removed.
#   Output: all_gape_with_mussel.csv, used as input for the next processing
#   step (valve gape and CV calculation)

library(tidyverse)


# 1. Load metadata

metadata <- read_delim(
  "data/raw/gape_metadata.csv",
  delim = ";",
  col_types = cols(.default = "c")
) %>%
  # Remove unused/empty window columns
  select(-c(`Gt1_1`, `Window6(min)`, `Gt1`, `Window1(min)`, `Gt4`, `Window4(min)`,
            `Gt5`, `Window5(min)`, `Gt2`, `Window2(min)`, `Gt3`, `Window3(min)`)) %>%
  mutate(
    Rep = as.numeric(Rep),
    Line_ID = as.numeric(Line_ID)
  )


View(metadata)

# 2. Set the CSV file path

path <- "data/raw/gape_raw"
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)

# 3. Define an auto CSV-reading function
# Detects whether the file is ";" or "," delimited, and cleans column
# names of BOM characters/whitespace

auto_read_csv <- function(file_path) {
  first_line <- readLines(file_path, n = 1)
  sep <- ifelse(grepl(";", first_line), ";", ",")
  df <- read_delim(file_path, delim = sep, show_col_types = FALSE)
  
  colnames(df) <- str_replace_all(colnames(df), "\ufeff", "") %>% str_trim()
  
  return(df)
}

# 4. Read all CSV files and reshape to long format

all_data <- map_df(files, function(f) {
  
  fname <- basename(f)
  
  # Extract Sensor_ID (A/B) and Rep number from the filename
  sensor <- str_extract(fname, "^(A|B)")
  rep <- str_extract(fname, "(?i)(?<=rep)\\d+") %>% as.numeric()
  
  print(paste("Reading file:", fname, "Sensor_ID:", sensor, "Rep:", rep))
  
  # Read CSV with automatic delimiter detection
  df <- auto_read_csv(f)
  
  # Standardise shared columns to consistent types across files
  if ("Day" %in% colnames(df)) {
    df <- df %>% mutate(Day = as.character(Day))
  } else {
    df$Day <- NA_character_
  }
  
  if ("Year" %in% colnames(df)) {
    df <- df %>% mutate(Year = as.numeric(Year))
  } else {
    df$Year <- NA_real_
  }
  
  if ("Month" %in% colnames(df)) {
    df <- df %>% mutate(Month = as.numeric(Month))
  } else {
    df$Month <- NA_real_
  }
  
  if ("Second" %in% colnames(df)) {
    df <- df %>% mutate(Second = as.numeric(Second))
  } else {
    df$Second <- NA_real_
  }
  
  if ("Minute" %in% colnames(df)) {
    df <- df %>% mutate(Minute = as.numeric(Minute))
  } else {
    df$Minute <- NA_real_
  }
  
  # Identify all channel columns (Ch1, Ch2, ... case-insensitive)
  ch_cols <- grep("(?i)^Ch\\d+", colnames(df), value = TRUE)
  print(paste("File:", fname, "Channel columns found:", paste(ch_cols, collapse = ", ")))
  
  if (length(ch_cols) == 0) {
    warning(paste("No channel columns found in file:", fname))
    return(NULL)
  }
  
  # Reshape to long format: one row per channel per timestamp
  df_long <- df %>%
    pivot_longer(
      cols = all_of(ch_cols),
      names_to = "Channel",
      values_to = "Gape"
    ) %>%
    mutate(
      Line_ID = as.numeric(str_remove(Channel, "(?i)Ch")),
      Sensor_ID = sensor,
      Rep = rep
    )
  
  # Match each channel to its corresponding mussel via metadata
  df_long %>%
    left_join(metadata, by = c("Sensor_ID", "Rep", "Line_ID"))
})

View(all_data)

# 5. Remove rows with no matching Mussel_ID (unused/empty channels)

all_data_clean <- all_data %>%
  filter(!is.na(Mussel_ID))

# 6. Export results

output_path <- "data/processed/all_gape_with_mussel.csv"
write_csv(all_data_clean, output_path)

View(all_data_clean)
cat("Processing complete. Rows with NA Mussel_ID removed. Result saved to:", output_path, "\n")
