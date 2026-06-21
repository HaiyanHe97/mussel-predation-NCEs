# Heart Rate Calculation - Raw PULSE Signal Processing
# Project: Effects of parasitic infection on predation risk responses in
#          blue mussels (Mytilus edulis)
# Description:
#   This script processes raw PULSE sensor recordings into heart rate
#   estimates using the heartbeatr package. For each batch folder:
#     1) Raw signal is read and split into 30-sec windows, sampled every
#        60 sec (window_width_secs = 30, window_shift_secs = 60)
#     2) Signal is optimised (interpolated and smoothed) for visual QC,
#        comparing raw vs optimised signal per channel
#     3) Heart rate is estimated via peak detection (pulse_heart) and
#        double-checked (pulse_doublecheck)
#     4) Heart rate estimates are summarised into 5-min windows using the
#        median (pulse_summarise default FUN), which is robust to
#        occasional abnormal readings caused by signal noise
#   Output: one heart rate CSV and two diagnostic plots (optimised vs raw
#   signal, heart rate pattern) per batch folder, saved to
#   data/processed/heart_rate_by_batch/
#   A separate "check" section at the end allows visual inspection of the
#   raw waveform for a specific individual/timestamp, used for spot-
#   checking peak detection accuracy.

library(heartbeatr)
library(dplyr)
library(stringr)
library(ggplot2)
library(fs)
library(tidyr)

# 1. Get the data file paths

files <- list.files(path = "data/raw/heart_rate_raw",
                    pattern = "\\.CSV$", recursive = TRUE, full.names = TRUE)
folders <- unique(path_dir(files))
pulse_multi <- TRUE

# 2. Process each batch folder: split, optimise, calculate heart rate

for (f in folders) {
  cat("Processing:", f, "\n")
  
  pulse_data_df_test <- pulse_read(paths = f, msg = T)
  
  # Split into time windows
  pulse_split_data <- pulse_split(pulse_data_df_test,
                                  window_width_secs = 30,
                                  window_shift_secs = 60)
  
  # Optimise (interpolate + smooth) for visual QC
  pulse_data_optimized <- pulse_optimize(
    pulse_data_split = pulse_split_data,
    interpolation_freq = 40,
    bandwidth = 0.75,
    raw_v_smoothed = FALSE,
    multi = TRUE
  )
  
  # Plot optimised data vs raw data (visual QC)
  pulse_plot <- bind_rows(
    pulse_split_data %>% mutate(type = "raw"),
    pulse_data_optimized %>% mutate(type = "opt")
  ) %>%
    unnest(data)
  
  pulse_long <- pulse_plot %>%
    pivot_longer(
      cols = starts_with("mussel_"),
      names_to = "mussel",
      values_to = "value"
    )
  
  p <- ggplot(pulse_long,
              aes(x = time,
                  y = value,
                  group = interaction(i, mussel, type))) +
    geom_line(
      data = subset(pulse_long, type == "raw"),
      color = "black",
      alpha = 0.5
    ) +
    geom_line(
      data = subset(pulse_long, type == "opt"),
      color = "red",
      linewidth = 0.8
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_blank()
    )
  
  ggsave(filename = paste0("data/processed/heart_rate_by_batch/", basename(f), "_optimized.jpg"),
         plot = p, width = 10, height = 4)
  
  # Calculate heart rates
  heart_rates <- pulse_heart(pulse_split_data, msg = FALSE)
  heart_rates <- pulse_doublecheck(heart_rates = heart_rates)
  
  # Summarise into 5-min windows (median, the pulse_summarise default)
  heart_rates_summary <- pulse_summarise(heart_rates, span_mins = 5)
  
  write.csv(heart_rates_summary,
            file = paste0("data/processed/heart_rate_by_batch/", basename(f), "_heartrate.csv"))
  
  p2 <- pulse_plot(heart_rates)
  ggsave(filename = paste0("data/processed/heart_rate_by_batch/", basename(f), "_hr_pattern.jpg"),
         plot = p2, width = 8, height = 8)
}

# 3. Spot-check: visually inspect peak detection for a specific batch

files <- list.files(path = "data/raw/heart_rate_raw",
                    pattern = "\\.CSV$", recursive = TRUE, full.names = TRUE)
folders <- unique(path_dir(files))
check <- folders[25]

pulse_data_df_check <- pulse_read(paths = check, msg = T)

pulse_split_data <- pulse_split(pulse_data_df_check,
                                window_width_secs = 30,
                                window_shift_secs = 60)

pulse_data_optimized <- pulse_optimize(
  pulse_data_split = pulse_split_data,
  interpolation_freq = 40,
  bandwidth = 0.75,
  raw_v_smoothed = FALSE,
  multi = TRUE
)

heart_rates <- pulse_heart(pulse_split_data, msg = FALSE)
heart_rates <- pulse_doublecheck(heart_rates = heart_rates)

head(heart_rates)
heart_rates$id
heart_rates$time

pulse_plot_raw(heart_rates, ID = "mussel_252", target_time = "2025-06-27 14:34:15.017 UTC")
