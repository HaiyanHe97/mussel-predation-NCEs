# =============================================================================
# Net Displacement Calculation - From Tracked Coordinates
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   This script calculates the net displacement (straight-line distance
#   between start and end position) for each individual mussel, based on
#   lens-distortion-corrected X/Y coordinate tracking data (one file per
#   arena replicate). Pixel distances are converted to mm using arena-
#   specific pixel-to-cm calibration values.
#   Output: net_distance_ALL.csv, used as input for 02_build_arena_model.R
# =============================================================================

library(dplyr)
library(stringr)
library(tidyr)

# 1. Function: calculate net displacement for one arena file

calc_net_distance <- function(file_path, px_per_cm) {
  
  total <- read.csv(
    file_path,
    header = TRUE,
    sep = ";",
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  
  cols <- colnames(total)
  
  # Identify X and Y coordinate columns (one pair per mussel)
  x_cols <- cols[grepl("^X_", cols)]
  y_cols <- cols[grepl("^Y_", cols)]
  
  # ---- Calculate net displacement (px) per mussel ----
  # Net displacement = straight-line distance between first and last
  # recorded position
  net_distance <- data.frame(
    Mussel = x_cols,
    NetDistance_px = NA
  )
  
  for (i in seq_along(x_cols)) {
    x_start <- total[[x_cols[i]]][1]
    y_start <- total[[y_cols[i]]][1]
    x_end   <- total[[x_cols[i]]][nrow(total)]
    y_end   <- total[[y_cols[i]]][nrow(total)]
    
    net_distance$NetDistance_px[i] <- sqrt(
      (x_end - x_start)^2 + (y_end - y_start)^2
    )
  }
  
  # ---- Parse mussel identity from column name (e.g. "X_Arena1_1_CM_3") ----
  df_name <- as.data.frame(
    str_split_fixed(net_distance$Mussel, "_", 5)
  )[, -1]
  
  colnames(df_name) <- c("Arena", "Rep", "Group", "ID")
  
  df_name$Arena <- as.numeric(sub("Arena", "", df_name$Arena))
  
  # ---- Merge displacement values ----
  df_name$NetDistance_px <- net_distance$NetDistance_px
  
  # ---- Convert px to mm using arena-specific calibration ----
  df_name$NetDistance_mm <- 10 * df_name$NetDistance_px / px_per_cm
  
  # Track source file (arena replicate identifier, e.g. "6A")
  df_name$File <- basename(file_path)
  
  return(df_name)
}

# 2. Apply function across all arena coordinate files

folder <- "data/raw/corrected_coordinates"

files <- list.files(
  folder,
  pattern = "Corrected.csv",
  full.names = TRUE
)

# ---- Pixel-to-cm calibration values, one per arena replicate file ----
px_info <- data.frame(
  File = basename(files),
  px_per_cm = c(
    26.5,   # 1A
    27.5,   # 2A
    26.67,  # 3A
    28.02,  # 4A
    26.33,  # 5A
    28.07   # 6A
  )
)

all_results <- lapply(files, function(f) {
  
  px_cm <- px_info$px_per_cm[px_info$File == basename(f)]
  
  calc_net_distance(
    file_path = f,
    px_per_cm = px_cm
  )
})

final_table <- bind_rows(all_results)

View(final_table)

# 3. Export results

write.csv(
  final_table,
  "data/processed/net_distance_ALL.csv",
  row.names = FALSE
)
