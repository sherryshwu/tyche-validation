#!/usr/bin/env Rscript
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
})

source("scripts/analysis/tree_functions.R")

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  cat("Usage: Rscript publication_plots.R summary_files simulation_names analysis_type rev_suffixes selected_models selected_configs plot_type output_dir [skip_convergence]\n")
  cat("plot_type options:\n")
  cat("  'main' - both primary and GC re-entry simulations, 1:1 configs, three metrics (height/RF/MRCA)\n")
  cat("  'supp_all_metrics' - both primary and GC re-entry simulations, 1:3 configs, three metrics\n") 
  cat("  'supp_tree_length' - both primary and GC re-entry simulations, all configs, tree length analysis\n")
  quit(status = 1)
}

summary_files <- strsplit(args[1], ",")[[1]]
simulation_names <- strsplit(args[2], ",")[[1]]
analysis_type <- args[3]
rev_suffixes <- strsplit(args[4], ",")[[1]]
selected_models <- strsplit(args[5], ",")[[1]]
selected_configs <- strsplit(args[6], ",")[[1]]
plot_type <- args[7]  # "main" or "supp_all_metrics" or "supp_tree_length"
output_dir <- args[8]
skip_convergence_filtering <- length(args) >= 9 && args[9] == "skip_convergence"

cat("=== Creating", toupper(plot_type), "Publication Plots ===\n")
cat("Summary files:", paste(summary_files, collapse=", "), "\n")
cat("Simulations:", paste(simulation_names, collapse=", "), "\n")
cat("Plot type:", plot_type, "\n")
cat("Rev suffixes:", paste(rev_suffixes, collapse=", "), "\n")
cat("Models:", paste(selected_models, collapse=", "), "\n")
cat("Configs:", paste(selected_configs, collapse=", "), "\n")
cat("Output directory:", output_dir, "\n")

if (skip_convergence_filtering) {
  cat("NOTE: Skipping convergence filtering - plotting all clones\n")
}

# Setup output directories
project_root <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
plots_dir <- file.path(output_dir, "plots")
pub_plots_dir <- file.path(plots_dir, "main_figures")
supp_plots_dir <- file.path(plots_dir, "supplementary_figures")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pub_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_plots_dir, showWarnings = FALSE)

# -------------------------- Load all data --------------------------
cat("Loading data from summary files...\n")

# Load combined summary
all_combined_data <- data.frame()
for (summary_file in summary_files) {
  if (file.exists(summary_file)) {
    temp_data <- read_csv(summary_file, show_col_types = FALSE)
    all_combined_data <- bind_rows(all_combined_data, temp_data)
    cat("Loaded", nrow(temp_data), "rows from", basename(summary_file), "\n")
  } else {
    cat("Warning: Summary file not found:", summary_file, "\n")
  }
}

if (nrow(all_combined_data) == 0) {
  cat("ERROR: No data loaded from summary files\n")
  quit(status = 1)
}

cat("Total combined data rows:", nrow(all_combined_data), "\n")

# -------------------------- Load convergence data --------------------------
conv_exclusions_by_config <- data.frame(config = character(), clone_id = character(), 
                                        simulation_name = character(), rev_suffix = character())

if (!skip_convergence_filtering) {
  cat("Loading convergence data...\n")
  
  # Build convergence file paths based on simulation names
  for (sim_name in simulation_names) {
    for (rev_suffix in rev_suffixes) {
      base_dir <- file.path(project_root, sim_name, "results", analysis_type, rev_suffix)      
      convergence_files <- c(
        file.path(base_dir, "tyche_models", "summary", "tyche_models_convergence_summary.csv"),
        file.path(base_dir, "competing_models", "summary", "competing_models_convergence_summary.csv")
      )

      for (conv_file in convergence_files) {
        if (file.exists(conv_file)) {
          temp_conv <- read_csv(conv_file, show_col_types = FALSE)
          temp_conv$unconverged_clone_ids <- as.character(temp_conv$unconverged_clone_ids)

          # Add tracking columns
          temp_conv$simulation_name <- sim_name
          temp_conv$rev_suffix <- rev_suffix

          # Process exclusions for this specific combination
          filtered_convergence_data <- temp_conv %>%
            mutate(model_short = get_model_short_name(template_id)) %>%
            filter(model_short %in% selected_models) %>%
            filter(config_name %in% selected_configs)

          if (nrow(filtered_convergence_data) > 0) {
            parse_clone_ids <- function(x) {
              if (is.na(x) || x == "") return(character(0))
              x <- gsub("[^0-9]", " ", x)
              parts <- unlist(strsplit(x, "\\s+"))
              parts <- parts[nzchar(parts)]
              unique(as.character(as.integer(parts)))
            }

            for (j in seq_len(nrow(filtered_convergence_data))) {
              ids <- parse_clone_ids(filtered_convergence_data$unconverged_clone_ids[j])
              if (length(ids) > 0) {
                temp_exclusions <- data.frame(
                  config = filtered_convergence_data$config_name[j], 
                  clone_id = ids,
                  simulation_name = sim_name,
                  rev_suffix = rev_suffix,
                  stringsAsFactors = FALSE
                )
                conv_exclusions_by_config <- rbind(conv_exclusions_by_config, temp_exclusions)
              }
            }
          }
          cat("Processed convergence data from:", basename(conv_file), "\n")
        }
      }
    }
  }

  conv_exclusions_by_config <- unique(conv_exclusions_by_config)
  cat("Total convergence exclusions:", nrow(conv_exclusions_by_config), "\n")
} else {
  cat("Skipping convergence filtering - plotting all clones...\n")
}

# Define model labels
model_labels <- c(
  "EO_Fixed" = "TyCHE\nfixed clocks",
  "EO_Est" = "TyCHE\nest clocks",
  "IS_Est" = "TyCHE (IS)\nest clocks",
  "MS_Fixed" = "TyCHE (MS)\nfixed clocks",
  "MS_Est" = "TyCHE (MS)\nest clocks",
  "SC_AR" = "SC",
  "UCLD_AR" = "UCLD"
)

# -------------------------- Process and filter data --------------------------
if (nrow(all_combined_data) > 0) {
  # Filter data based on selections
  filtered_data_initial <- all_combined_data %>%
    filter(
      simulation_name %in% simulation_names,
      rev_suffix %in% rev_suffixes
    )

  # Make join keys comparable and filter out unconverged clones
  filtered_data_initial <- filtered_data_initial %>%
    mutate(clone_id = as.character(clone_id))
  conv_exclusions_by_config <- conv_exclusions_by_config %>%
    mutate(clone_id = as.character(clone_id))

  filtered_data <- filtered_data_initial %>%
    anti_join(conv_exclusions_by_config, 
              by = c("config", "clone_id", "simulation_name", "rev_suffix"))

  dropped_n <- nrow(filtered_data_initial) - nrow(filtered_data)
  cat("Filtered out", dropped_n, "rows from unconverged clones\n")

  # Process the data for plotting
  summary_data <- filtered_data %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(
      # Calculate metrics
      tree_height_prop_error = (beast_tree_height - true_tree_height) / true_tree_height,
      tree_length_prop_error = (beast_tree_length - true_tree_length) / true_tree_length,
      mrca_cell_type_accuracy = 100 * mrca_cell_type_accuracy,
      clone_id = factor(clone_id, levels = as.character(1:20)),
      model_short = get_model_short_name(template_name)
    ) %>%
    filter(model_short %in% selected_models) %>%
    mutate(
      model_short = factor(model_short, levels = selected_models),
      model_display = factor(model_labels[as.character(model_short)],
                             levels = model_labels[selected_models]),
      # Create facet labels
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
      )),
      # Create basic simulation-specific labels (for main and supp_all_metrics)
      facet_label_with_sim = if (plot_type %in% c("main", "supp_all_metrics")) {
        case_when(
          simulation_name == "gc_reentry_hunter" & grepl("_sel$", config) ~ "Selective Evolution (GC re-entry)",
          simulation_name == "gc_reentry_hunter" & grepl("_neu$", config) ~ "Uniform Neutral Evolution (GC re-entry)",
          simulation_name == "tltt_08_20" & grepl("_sel$", config) ~ "Selective Evolution",
          simulation_name == "tltt_08_20" & grepl("_neu$", config) ~ "Uniform Neutral Evolution",
          TRUE ~ as.character(facet_label)
        )
      } else {
        as.character(facet_label)
      },

      # Create detailed labels with config info (for supp_tree_length)
      facet_label_detailed = if (plot_type == "supp_tree_length") {
        case_when(
          simulation_name == "gc_reentry_hunter" & config == "config_ratio_1to1_sel" ~ "Selective Evolution (GC re-entry)\n1:1 GC:Other",
          simulation_name == "gc_reentry_hunter" & config == "config_ratio_1to1_neu" ~ "Uniform Neutral Evolution (GC re-entry)\n1:1 GC:Other",
          simulation_name == "gc_reentry_hunter" & config == "config_ratio_1to3_sel" ~ "Selective Evolution (GC re-entry)\n1:3 GC:Other",
          simulation_name == "gc_reentry_hunter" & config == "config_ratio_1to3_neu" ~ "Uniform Neutral Evolution (GC re-entry)\n1:3 GC:Other",
          simulation_name == "tltt_08_20" & config == "config_ratio_1to1_sel" ~ "Selective Evolution\n1:1 GC:Other",
          simulation_name == "tltt_08_20" & config == "config_ratio_1to1_neu" ~ "Uniform Neutral Evolution\n1:1 GC:Other",
          simulation_name == "tltt_08_20" & config == "config_ratio_1to3_sel" ~ "Selective Evolution\n1:3 GC:Other",
          simulation_name == "tltt_08_20" & config == "config_ratio_1to3_neu" ~ "Uniform Neutral Evolution\n1:3 GC:Other",
          TRUE ~ as.character(facet_label)
        )
      } else {
        as.character(facet_label)
      },

      # Set factor levels for each
      facet_label_with_sim = factor(facet_label_with_sim, levels = c(
        "Selective Evolution", "Uniform Neutral Evolution",
        "Selective Evolution (GC re-entry)", "Uniform Neutral Evolution (GC re-entry)"
      )),

      facet_label_detailed = factor(facet_label_detailed, levels = c(
        "Selective Evolution\n1:1 GC:Other", 
        "Uniform Neutral Evolution\n1:1 GC:Other",
        "Selective Evolution\n1:3 GC:Other", 
        "Uniform Neutral Evolution\n1:3 GC:Other",
        "Selective Evolution (GC re-entry)\n1:1 GC:Other", 
        "Uniform Neutral Evolution (GC re-entry)\n1:1 GC:Other",
        "Selective Evolution (GC re-entry)\n1:3 GC:Other", 
        "Uniform Neutral Evolution (GC re-entry)\n1:3 GC:Other"
      ))
    ) %>%
    filter(config %in% selected_configs)  # Apply config filter
  
  cat("After filtering:", nrow(summary_data), "data rows for plotting\n")
}

# -------------------------- Main publication figures --------------------------
plots <- list()
# Determine output directory and facet column
if (plot_type == "main") {
  output_figure_dir <- pub_plots_dir
  facet_column <- "facet_label_with_sim"
} else if (plot_type == "supp_all_metrics") {
  output_figure_dir <- supp_plots_dir
  facet_column <- "facet_label_with_sim"
} else if (plot_type == "supp_tree_length") {
  output_figure_dir <- supp_plots_dir
  facet_column <- "facet_label_detailed"
}

if (plot_type %in% c("main", "supp_all_metrics")) {
  cat("Creating three-metric publication figures...\n")
  cat("Target for main: Two simulations x Two 1:1 configs = Four panels\n")
  cat("Target for supp_all_metrics: Two simulations x Two 1:3 configs = Four panels\n")
  cat("Using facet column:", facet_column, "\n")

  # Create individual plots
  if (length(summary_data$tree_height_prop_error) > 0) {
    plots$height <- save_metric_plot(
      metric_col = "tree_height_prop_error",
      y_label = "Tree Height\nProportional Error",
      show_x_labels = FALSE,
      add_reference_line = TRUE,
      reference_value = 0,
      facet_col = facet_column
    ) + theme(
      text = element_text(size = 8),
      strip.text = element_text(size = 9)  # Keep strip text for top row
    )
  }

  if (length(summary_data$rf_distance) > 0) {
    plots$rf_distance <- save_metric_plot(
      metric_col = "rf_distance", 
      y_label = "RF Distance",
      show_x_labels = FALSE,
      add_reference_line = TRUE,
      reference_value = 0,
      facet_col = facet_column
    ) + theme(
      text = element_text(size = 8),
      strip.text = element_blank()  # Remove strip text for middle row
    )
  }

  if (length(summary_data$mrca_cell_type_accuracy) > 0) {
    plots$mrca <- save_metric_plot(
      metric_col = "mrca_cell_type_accuracy",
      y_label = "Ancestral Cell\nType Accuracy (%)",
      show_x_labels = TRUE,
      add_reference_line = TRUE,
      reference_value = 100,
      facet_col = facet_column
    ) + theme(
      text = element_text(size = 8),
      strip.text = element_blank()  # Remove strip text for bottom row
    )
  }

  # Save individual plots
  if (!is.null(plots$height)) {
    ggsave(file.path(output_figure_dir, paste0(plot_type, "_tree_height.pdf")),
           plots$height, width = 12, height = 8)
  }

  if (!is.null(plots$rf_distance)) {
    ggsave(file.path(output_figure_dir, paste0(plot_type, "_rf_distance.pdf")),
           plots$rf_distance, width = 12, height = 8)
  }

  if (!is.null(plots$mrca)) {
    ggsave(file.path(output_figure_dir, paste0(plot_type, "_mrca_accuracy.pdf")),
           plots$mrca, width = 12, height = 8)
  }

  # Create combined plot
  if (length(plots) >= 3) {
    combined_plot <- plots$height / plots$rf_distance / plots$mrca +
      plot_layout(heights = c(1, 1, 1)) &
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(1, 2, 1, 2) # t, r, b, l
      )

    ggsave(file.path(output_figure_dir, paste0(plot_type, "_combined_metrics.pdf")),
           combined_plot, width = 7, height = 3.5)

    cat("✓ Main publication plots saved\n")
  }

} else if (plot_type == "supp_tree_length") {
  cat("Creating supplementary figures for tree length...\n")
  cat("Target: Both simulations x Four configs = Eight panels\n")
  cat("Using facet column:", facet_column, "\n")

  # Tree length proportional error analysis
  if (length(summary_data$tree_length_prop_error) > 0) {
    tree_length_plot <- save_metric_plot(
      metric_col = "tree_length_prop_error",
      y_label = "Tree Length\nProportional Error",
      show_x_labels = TRUE,
      add_reference_line = TRUE,
      reference_value = 0,
      facet_col = facet_column,
      nrow = 2,
      ncol = 4
    ) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(1, 2, 1, 2) # t, r, b, l
      )

    ggsave(file.path(output_figure_dir, paste0(plot_type, "_prop_error.pdf")),
           tree_length_plot, width = 10, height = 5)
  }
  cat("✓ Supplementary plots saved\n")
}

# -------------------------- Final summary --------------------------
cat("\n=== Publication Plots Complete ===\n")
cat("Plot type:", plot_type, "\n")
cat("Output directory:", pub_plots_dir, "\n")

# List all created files
pdf_files <- list.files(pub_plots_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = FALSE)

cat("\nTotal files created:", length(pdf_files), "\n")

cat("\n=== Key Publication Files ===\n")
cat("Files created:\n")
for (file in pdf_files) {
  cat("  -", file, "\n")
}

cat("=== Done ===\n")