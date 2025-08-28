#!/usr/bin/env Rscript
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
})

source("scripts/tree_functions.R")

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript publication_plots.R tree_analysis_dir simulation_name analysis_type rev_suffix\n")
  quit(status = 1)
}

tree_analysis_dir <- args[1]
simulation_name <- args[2]
analysis_type <- args[3]
rev_suffix <- args[4]
selected_models <- strsplit(args[5], ",")[[1]]

cat("=== Creating Publication Plots ===\n")
cat("Analysis directory:", tree_analysis_dir, "\n")

# Setup output directories
plots_dir <- file.path(tree_analysis_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
pub_plots_dir <- file.path(plots_dir, "publication")
dir.create(pub_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(pub_plots_dir, "supplementary"), showWarnings = FALSE)

# -------------------------- Load all data --------------------------
cat("Loading data...\n")

# Load combined summary
summary_file <- file.path(tree_analysis_dir, "all_results_summary.csv")
if (!file.exists(summary_file)) {
  cat("ERROR: Combined summary not found:", summary_file, "\n")
  quit(status = 1)
}

all_data <- read_csv(summary_file, show_col_types = FALSE)
cat("Loaded", nrow(all_data), "data rows\n")

# -------------------------- Load convergence data --------------------------
cat("Loading convergence summaries...\n")

# Load convergence summaries for different model types
convergence_files <- c(
  file.path(dirname(tree_analysis_dir), "tyche_models", "summary", "tyche_models_convergence_summary.csv"),
  file.path(dirname(tree_analysis_dir), "competing_models", "summary", "competing_models_convergence_summary.csv")
)

convergence_data <- data.frame()
for (conv_file in convergence_files) {
  if (file.exists(conv_file)) {
    temp_conv <- read_csv(conv_file, show_col_types = FALSE)
    temp_conv$unconverged_clone_ids <- as.character(temp_conv$unconverged_clone_ids)
    convergence_data <- bind_rows(convergence_data, temp_conv)
    cat("Loaded convergence data from:", basename(conv_file), "\n")
  } else {
    cat("Warning: Convergence file not found:", conv_file, "\n")
  }
}

if (nrow(convergence_data) == 0) {
  cat("ERROR: No convergence data found. Cannot filter unconverged runs.\n")
  quit(status = 1)
}

cat("Building clone-level exclusions from convergence summaries...\n")

parse_clone_ids <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  x <- gsub("[^0-9]", " ", x)                  # replace non-digits with spaces
  parts <- unlist(strsplit(x, "\\s+"))
  parts <- parts[nzchar(parts)]
  unique(as.character(as.integer(parts)))      # normalize like "01" -> "1"
}

# Expand into (config, template, clone_id) rows using a simple loop
conv_rows <- list()
if (nrow(convergence_data) > 0) {
  for (i in seq_len(nrow(convergence_data))) {
    ids <- parse_clone_ids(as.character(convergence_data$unconverged_clone_ids[i]))
    if (length(ids) == 0) next
    conv_rows[[length(conv_rows) + 1]] <- data.frame(
      config        = as.character(convergence_data$config_name[i]),
      template_name = as.character(convergence_data$template_id[i]),
      clone_id      = ids,
      stringsAsFactors = FALSE
    )
  }
}

conv_exclusions <- if (length(conv_rows) == 0) {
  data.frame(config = character(), template_name = character(), clone_id = character(), stringsAsFactors = FALSE)
} else {
  unique(do.call(rbind, conv_rows))
}

conv_summary <- convergence_data %>%
  group_by(config_name, template_id, model_type) %>%
  summarise(total_clones = first(total_clones),
            converged_clones = first(converged_clones),
            .groups = "drop")
cat("\nConvergence Summary (per config-template-model):\n")
print(conv_summary)

model_labels <- c(
  "EO_Fixed" = "TyCHE\nfixed clocks",
  "EO_Est" = "TyCHE\nest clocks",
  "IS_Est" = "TyCHE (IS)\nest clocks",
  "MS_Fixed" = "TyCHE (MS)\nfixed clocks",
  "MS_Est" = "TyCHE (MS)\nest clocks",
  "SC_AR" = "SC",
  "UCLD_AR" = "UCLD"
)

if (nrow(all_data) > 0) {
  # Make join keys comparable
  all_data <- all_data %>% mutate(clone_id = as.character(clone_id))
  conv_exclusions <- conv_exclusions %>% mutate(clone_id = as.character(clone_id))

  # Drop only the unconverged clones via anti_join
  filtered_data <- all_data %>%
    anti_join(conv_exclusions, by = c("config", "template_name", "clone_id"))

  dropped_n <- nrow(all_data) - nrow(filtered_data)
  if (dropped_n > 0) {
    cat("Filtered out", dropped_n, "rows from unconverged clone IDs\n")
  } else {
    cat("No rows matched unconverged clone IDs; nothing dropped\n")
  }

  summary_data <- filtered_data %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    # Calculate tree height proportional error
    mutate(prop_error = (beast_tree_height - true_tree_height) / true_tree_height) %>%
    # Extract abbreviations for model names
    mutate(model_short = get_model_short_name(template_name)) %>%
    filter(model_short %in% selected_models) %>%
    mutate(
      mrca_cell_type_accuracy = 100 * mrca_cell_type_accuracy,
      clone_id = factor(clone_id, levels = as.character(1:20)),
      model_short = factor(model_short, levels = selected_models),
      model_display = factor(model_labels[as.character(model_short)], levels = model_labels[selected_models])
    ) %>%
    # Extract ratio info for facet grid
    mutate(
      facet_label = case_when(
        grepl("ratio_1to3", config) & grepl("_sel$", config) ~ "Selective Evolution\n1:3 GC:Other",
        grepl("ratio_1to3", config) & grepl("_neu$", config) ~ "Uniform Neutral Evolution\n1:3 GC:Other",
        grepl("ratio_1to1", config) & grepl("_sel$", config) ~ "Selective Evolution\n1:1 GC:Other",
        grepl("ratio_1to1", config) & grepl("_neu$", config) ~ "Uniform Neutral Evolution\n1:1 GC:Other",
        TRUE ~ "other"
      ),
      facet_label = factor(facet_label, levels = c(
        "Selective Evolution\n1:3 GC:Other",
        "Uniform Neutral Evolution\n1:3 GC:Other",
        "Selective Evolution\n1:1 GC:Other",
        "Uniform Neutral Evolution\n1:1 GC:Other"
      ))
    )
  cat("After filtering for convergence:", nrow(summary_data), "data rows remaining\n")
}

# -------------------------- Main publication figures --------------------------
cat("Creating main figures...\n")

plots <- list()
if (length(summary_data$prop_error) > 0) {
  # Figure 2C: Tree height comparison
  cat("Creating Figure 2C: Tree Heights\n")
  tree_height_plot <- save_metric_plot(
    metric_col = "beast_tree_height",
    y_label = "Tree Height",
    title = "Estimated vs True Tree Heights",
    show_x_labels = TRUE,
    add_reference_line = TRUE,
    reference_value = median(summary_data$true_tree_height)
  )

  ggsave(file.path(pub_plots_dir, "tree_heights.pdf"),
         tree_height_plot, width = 12, height = 8)

  plots$height <- save_metric_plot(
    metric_col = "prop_error",
    y_label = "Tree Height\n Proportional Error",
    show_x_labels = FALSE,
    add_reference_line = TRUE,
    reference_value = 0
  )
}

if (length(summary_data$rf_distance) > 0) {
  # Figure 2D: RF Distance
  cat("Creating Figure 2D: RF Distance\n")
  plots$rf_distance <- save_metric_plot(
    metric_col = "rf_distance",
    y_label = "RF Distance",
    show_x_labels = FALSE,
    add_reference_line = TRUE,
    reference_value = 0
  )
  ggsave(file.path(pub_plots_dir, "rf_distance.pdf"),
         plots$rf_distance, width = 12, height = 8)
}

if (length(summary_data$mrca_cell_type_accuracy) > 0) {
  # Figure 2E: Ancestral Cell Type Accuracy
  cat("Creating Figure 2E: Ancestral Cell Type Accuracy\n")
  plots$mrca <- save_metric_plot(
    metric_col = "mrca_cell_type_accuracy",
    y_label = "Ancestral Cell\n Type Accuracy (%)",
    show_x_labels = TRUE,
    add_reference_line = TRUE,
    reference_value = 100
  )
  ggsave(file.path(pub_plots_dir, "mrca_cell_type_accuracy.pdf"),
         plots$mrca, width = 12, height = 8)
}

# Apply compact theme and remove strip text for middle and bottom rows
if (!is.null(plots$height)) {
  plots$height <- plots$height +
    theme(
      text = element_text(size = 8),
      strip.text = element_text(size = 9)  # Keep strip text for top row
    )
}

if (!is.null(plots$rf_distance)) {
  plots$rf_distance <- plots$rf_distance + theme(
    text = element_text(size = 8),
    strip.text = element_blank()  # Remove strip text for middle row
  )
}

if (!is.null(plots$mrca)) {
  plots$mrca <- plots$mrca + theme(
    text = element_text(size = 8),
    strip.text = element_blank()  # Remove strip text for bottom row
  )
}

if (!is.null(plots$height) && !is.null(plots$rf_distance) && !is.null(plots$mrca)) {
  combined_plot <- plots$height / plots$rf_distance / plots$mrca +
    plot_layout(heights = c(1, 1, 1)) &
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      plot.margin = margin(1, 2, 1, 2) # t, r, b, l
    )

  # Save using save_compact_plot with custom dimensions for 3x4 layout
  combined_path <- file.path(pub_plots_dir, "combined_metrics.pdf")
  ggsave(combined_path, combined_plot, width = 7, height = 3.5)

  cat("✓ Combined 3x4 faceted plot saved to:", combined_path, "\n")
} else {
  cat("Could not create combined plot - missing data for one or more metrics\n")
}

# -------------------------- Final summary --------------------------
cat("\n=== Publication Plots Complete ===\n")
cat("Output directory:", pub_plots_dir, "\n")

# List all created files
pdf_files <- list.files(pub_plots_dir, pattern = "\\.pdf$", recursive = TRUE)

cat("\nTotal files created:", length(pdf_files), "\n")

cat("\n=== Key Publication Files ===\n")
cat("Main figures (full size):\n")
cat("  - tree_heights.pdf\n")
cat("  - rf_distance.pdf\n")
cat("  - mrca_accuracy.pdf\n")

cat("\nCompact combined figure:\n")
cat("  - combined_metrics.pdf\n")

# cat("\nSpecialized analyses:\n")

# cat("\nSupplementary figures in: supplementary/\n")
# cat("=========================================\n")