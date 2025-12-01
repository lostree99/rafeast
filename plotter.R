#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Data Loading & Filtering
# -----------------------------------------------------------------------------
result_dir <- "results_final"
files <- list.files(result_dir, pattern = "*.csv", full.names = TRUE)

if (length(files) == 0) stop("No CSV files found!")

cat("Reading", length(files), "files...\n")

raw_data <- files %>%
  map_dfr(read_csv, show_col_types = FALSE)

# CRITICAL FIX: Truncate at N=16,000
# We exclude N > 16000 because Standard FEAST fails to converge (Error > 1.0),
# making the speedup comparison mathematically invalid beyond this point.
clean_data <- raw_data %>%
  filter(n <= 16000) %>%
  filter(!is.na(time_total))

cat("Data filtered to N <= 16,000. Proceeding with analysis.\n")

# -----------------------------------------------------------------------------
# 2. Pre-processing
# -----------------------------------------------------------------------------
summary_data <- clean_data %>%
  group_by(n) %>%
  mutate(
    baseline_time = mean(time_total[method_type == "Standard"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    speedup_instance = baseline_time / time_total,
    plot_label = case_when(
      method_type == "Standard" ~ "Standard FEAST (Nc=8)",
      method_type == "RA-FEAST" ~ paste0("RA-FEAST (Nc=", Nc, ", It=", max_iter, ")")
    )
  )

# Summarize for plotting
plot_df <- summary_data %>%
  group_by(n, plot_label, method_type, Nc, max_iter) %>%
  summarise(
    time_mean = mean(time_total, na.rm = TRUE),
    time_se   = sd(time_total, na.rm = TRUE) / sqrt(n()),
    speedup_mean = mean(speedup_instance, na.rm = TRUE),
    speedup_se   = sd(speedup_instance, na.rm = TRUE) / sqrt(n()),
    error_mean   = mean(error_vs_truth, na.rm = TRUE),
    p1_mean      = mean(time_phase1, na.rm = TRUE),
    p2_mean      = mean(time_phase2, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Plot 1: Speedup (The "Money Plot")
# -----------------------------------------------------------------------------
p1 <- ggplot(plot_df %>% filter(method_type == "RA-FEAST"), 
             aes(x = n, y = speedup_mean, color = plot_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  annotate("text", x = min(plot_df$n), y = 1.2, label = "Standard FEAST Baseline", color = "black", hjust = 0, vjust = 0) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = speedup_mean - speedup_se, ymax = speedup_mean + speedup_se), width = 0.1) +
  scale_x_continuous(trans = "log10", breaks = unique(plot_df$n)) +
  labs(title = "A. Computational Speedup", subtitle = "Relative to Standard FEAST", y = "Speedup Factor (X times)") +
  theme_bw() + theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Plot 2: Wall-Clock Time
# -----------------------------------------------------------------------------
p2 <- ggplot(plot_df, aes(x = n, y = time_mean, color = plot_label)) +
  geom_line(aes(linetype = method_type), linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(trans = "log10", breaks = unique(plot_df$n)) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c(scales::hue_pal()(length(unique(plot_df$plot_label))-1), "black")) +
  labs(title = "B. Wall-Clock Time Scaling", y = "Time (seconds, log scale)") +
  theme_bw() + theme(legend.position = "right")

# -----------------------------------------------------------------------------
# Plot 3: Overhead Analysis (At new Max N = 16,000)
# -----------------------------------------------------------------------------
max_n <- max(plot_df$n) # This will now be 16000
breakdown_df <- plot_df %>%
  filter(n == max_n) %>%
  select(plot_label, p1_mean, p2_mean) %>%
  pivot_longer(cols = c(p1_mean, p2_mean), names_to = "Phase", values_to = "Time") %>%
  mutate(
    Phase = ifelse(Phase == "p1_mean", "Phase 1: Warmstart", "Phase 2: Solver"),
    plot_label = reorder(plot_label, Time) 
  )

p3 <- ggplot(breakdown_df, aes(x = plot_label, y = Time, fill = Phase)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Phase 1: Warmstart" = "#D55E00", "Phase 2: Solver" = "#0072B2")) +
  labs(title = paste("C. Time Breakdown (N =", max_n, ")"), x = NULL, y = "Time (seconds)") +
  theme_bw() + theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Plot 4: Accuracy (Now Clean)
# -----------------------------------------------------------------------------
p4 <- ggplot(plot_df, aes(x = n, y = error_mean, color = plot_label)) +
  geom_line(aes(linetype = method_type), linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(trans = "log10", breaks = unique(plot_df$n)) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c(scales::hue_pal()(length(unique(plot_df$plot_label))-1), "black")) +
  labs(title = "D. Accuracy Validation", y = "Max Error vs Truth") +
  theme_bw() + theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Combine and Save
# -----------------------------------------------------------------------------
final_plot <- (p1 + p2) / (p3 + p4) + 
  plot_annotation(
    title = "RA-FEAST Performance Analysis",
    subtitle = "Standard FEAST (Black Line) vs. RA-FEAST Configurations (Up to N=16k)",
    caption = "Error bars represent Standard Error (SE)"
  )

ggsave("figure1_clean.png", final_plot, width = 14, height = 10, dpi = 300)
cat("Saved figure1_clean.png\n")