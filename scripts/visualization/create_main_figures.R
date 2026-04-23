#!/usr/bin/env Rscript
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
})

source("scripts/utils/phylo_utilities.R")

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  cat("Usage: Rscript create_main_figures.R summary_files simulation_names analysis_type rev_suffixes selected_models selected_configs plot_type output_dir [skip_convergence]\n")
  cat("plot_type options:\n")
  cat("  'main' - both primary and GC re-entry simulations, 1:1 sel config, four metrics (height/length/RF/MRCA)\n")
  cat("  'supp' - both primary and GC re-entry simulations, 1:1 neu config, four metrics (height/length/RF/MRCA)\n")
  quit(status = 1)
}

# -------------------------- Setup directories --------------------------
summary_files <- strsplit(args[1], ",")[[1]]
simulation_names <- strsplit(args[2], ",")[[1]]
analysis_type <- args[3]
rev_suffixes <- strsplit(args[4], ",")[[1]]
selected_models <- strsplit(args[5], ",")[[1]]
selected_configs <- strsplit(args[6], ",")[[1]]
plot_type <- args[7]  # "main" or "supp"
output_dir <- args[8]
skip_convergence_filtering <- length(args) >= 9 && args[9] == "skip_convergence"
exclude_clone_list <- if (length(args) >= 10) strsplit(args[10], ",")[[1]] else character(0)

cat("=== Creating", toupper(plot_type), "Publication Plots ===\n")
cat("Summary files:", paste(summary_files, collapse = ", "), "\n")
cat("Simulations:", paste(simulation_names, collapse = ", "), "\n")
cat("Plot type:", plot_type, "\n")
cat("Rev suffixes:", paste(rev_suffixes, collapse = ", "), "\n")
cat("Models:", paste(selected_models, collapse = ", "), "\n")
cat("Configs:", paste(selected_configs, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")

if (skip_convergence_filtering) {
  cat("NOTE: Skipping convergence filtering - plotting all clones\n")
}

# Setup output directories
project_root <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
plots_dir <- file.path(output_dir, "plots")
pub_plots_dir <- file.path(plots_dir, "main_figures")
supp_plots_dir <- file.path(plots_dir, "supplementary_figures")
plot_data_dir <- file.path(output_dir, "plot_data")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pub_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_plots_dir, showWarnings = FALSE)
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)

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

find_unconverged_clones_by_ess <- function(base_dir, config_name, selected_models, ess_threshold = 100) {
  all_exclusions <- data.frame(
    config = character(),
    clone_id = character(),
    model_abbr = character(),
    stringsAsFactors = FALSE
  )

  # Search only within the specific config directory
  config_dir <- file.path(base_dir, config_name)
  tree_files <- list.files(config_dir, pattern = "_tree_with_trait\\.tree$", recursive = TRUE, full.names = TRUE)

  if (length(tree_files) == 0) {
    cat("No tree files found in", basename(config_dir), "\n")
    return(all_exclusions)
  }

  cat("Found", length(tree_files), "tree files in", basename(config_dir), "\n")

  for (tree_file in tree_files) {
    filename <- basename(tree_file)

    # Extract model name from filename
    file_name_parts <- unlist(strsplit(filename, "_"))
    full_model_name <- paste(file_name_parts[1], file_name_parts[2], sep = "_")
    model_abbr <- abbreviate_model(full_model_name)

    # Process the tree if model is selected
    if (model_abbr %in% selected_models) {
      # Extract clone ID
      clone_id <- file_name_parts[which(file_name_parts == "tree") - 1]

      # Read ESS value
      ess_value <- read_ess_from_log(tree_file)

      if (!is.na(ess_value)) {
        if (ess_value <= ess_threshold) {
          # Add to exclusions
          exclusion <- data.frame(
            config = config_name,
            clone_id = clone_id,
            model_abbr = model_abbr,
            stringsAsFactors = FALSE
          )
          all_exclusions <- rbind(all_exclusions, exclusion)
          cat("Clone", clone_id, "for model", model_abbr, "ESS:", ess_value, "(unconverged)\n")
        }
      }
    }
  }
  return(all_exclusions)
}

# -------------------------- Load convergence data --------------------------
conv_exclusions_by_config <- data.frame(config = character(), clone_id = character(),
                                        simulation_name = character(), rev_suffix = character())

if (!skip_convergence_filtering) {
  cat("Loading convergence data...\n")

  # Build convergence file paths based on simulation names
  for (sim_name in simulation_names) {
    # Define simulation to rev_suffix mapping
    if (grepl("tltt", sim_name)) {
      rev_suffix <- "irrev"
    } else if (grepl("gc_reentry", sim_name)) {
      rev_suffix <- "rev"
    } else {
      cat("Warning: Could not determine rev_suffix for simulation:", sim_name, "- skipping convergence lookup\n")
      next
    }

    cat("Detected", rev_suffix, "model for simulation:", sim_name, "\n")

    base_dir <- file.path(project_root, sim_name, "results", analysis_type, rev_suffix)
    convergence_files <- c(
      file.path(base_dir, "tyche_models", "summary", "tyche_models_convergence_summary.csv"),
      file.path(base_dir, "competing_models", "summary", "competing_models_convergence_summary.csv")
    )

    # Check if convergence summary files exist
    csv_files_exist <- file.exists(convergence_files)

    if (any(csv_files_exist)) {
      cat("Using CSV convergence files for", sim_name, rev_suffix, "\n")

      for (conv_file in convergence_files) {
        if (file.exists(conv_file)) {
          temp_conv <- read_csv(conv_file, show_col_types = FALSE)
          temp_conv$unconverged_clone_ids <- as.character(temp_conv$unconverged_clone_ids)

          # Add tracking columns
          temp_conv$simulation_name <- sim_name
          temp_conv$rev_suffix <- rev_suffix

          # Process exclusions for this specific combination
          filtered_convergence_data <- temp_conv %>%
            mutate(model_abbr = abbreviate_model(template_id)) %>%
            filter(model_abbr %in% selected_models) %>%
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
    } else {
      cat("CSV convergence files missing for", sim_name, rev_suffix, "- using ESS method\n")

      # Use ESS method to find unconverged clones
      for (config_name in selected_configs) {
        for (model_abbr in selected_models) {
          # Look in both tyche_models and competing_models directories
          search_dirs <- c(
            file.path(base_dir, "tyche_models", "beast_raw_output"),
            file.path(base_dir, "competing_models", "beast_raw_output")
          )

          for (search_dir in search_dirs) {
            if (dir.exists(search_dir)) {
              exclusions_from_ess <- find_unconverged_clones_by_ess(search_dir, config_name, selected_models)

              if (nrow(exclusions_from_ess) > 0) {
                # Add simulation_name and rev_suffix columns
                exclusions_from_ess$simulation_name <- sim_name
                exclusions_from_ess$rev_suffix <- rev_suffix

                # Reorder columns to match the expected structure
                exclusions_from_ess <- exclusions_from_ess[, c("config", "clone_id", "simulation_name", "rev_suffix")]

                conv_exclusions_by_config <- rbind(conv_exclusions_by_config, exclusions_from_ess)
                cat("Found", nrow(exclusions_from_ess), "unconverged clones for", config_name, "using ESS method\n")

              }
            }
          }
        }
      }
    }
  }

  conv_exclusions_by_config <- unique(conv_exclusions_by_config)
  cat("Total convergence exclusions:", nrow(conv_exclusions_by_config), "\n")
} else {
  cat("Skipping convergence filtering - plotting all clones...\n")
}

# -------------------------- Process and filter data --------------------------
# Define model labels
model_labels <- c(
  "EO_Est" = "TyCHE",
  "SC_AR" = "SC",
  "UCLD_AR" = "UCLD"
)

if (nrow(all_combined_data) > 0) {
  # Filter data based on selections
  filtered_data_initial <- all_combined_data %>%
    filter(simulation_name %in% simulation_names, rev_suffix %in% rev_suffixes) %>%

    # Make join keys comparable and filter out unconverged clones
    mutate(
      clone_id = as.character(clone_id),
      model_abbr = abbreviate_model(template_name)
    ) %>%
    filter(model_abbr %in% selected_models)

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
      model_abbr = factor(model_abbr, levels = selected_models),
      model_display = factor(model_labels[as.character(model_abbr)],
                             levels = model_labels[selected_models]),
      # Create facet labels
      # Create basic simulation-specific labels (for main and supp)
      facet_label_with_sim = case_when(
        grepl("gc_reentry", simulation_name) & grepl("_sel$", config) ~ "Selective Evolution\n(GC re-entry)",
        grepl("gc_reentry", simulation_name) & grepl("_neu$", config) ~ "Uniform Neutral Evolution\n(GC re-entry)",
        grepl("tltt", simulation_name) & grepl("_sel$", config) ~ "Selective Evolution",
        grepl("tltt", simulation_name) & grepl("_neu$", config) ~ "Uniform Neutral Evolution",
        TRUE ~ as.character(paste0(simulation_name, config))
      ),

      # Set factor levels for each
      facet_label_with_sim = factor(facet_label_with_sim, levels = c(
        "Selective Evolution", "Uniform Neutral Evolution",
        "Selective Evolution\n(GC re-entry)", "Uniform Neutral Evolution\n(GC re-entry)"
      ))
    ) %>%
    filter(config %in% selected_configs)

  # Exclude specific clones if provided
  if (length(exclude_clone_list) > 0) {
    cat("Excluding specific clones:", paste(exclude_clone_list, collapse = ", "), "\n")
    summary_data <- summary_data %>% filter(!(clone_id %in% exclude_clone_list))
  }
  write_csv(summary_data, file.path(plot_data_dir, paste0("summary_data_", plot_type, ".csv")))
  cat("After filtering:", nrow(summary_data), "data rows for plotting\n")
}

# -------------------------- Main publication figures --------------------------
plots <- list()
facet_column <- "facet_label_with_sim"
# Determine output directory and facet column
if (plot_type == "main") {
  output_figure_dir <- pub_plots_dir
} else if (plot_type == "supp") {
  output_figure_dir <- supp_plots_dir
} else {
  cat("ERROR: Unknown plot type:", plot_type, "\n")
  quit(status = 1)
}

# -------------------------- Plotting logic --------------------------
create_main_panels <- function(data_subset, is_primary = TRUE) {
  # Helper to decide if we show the label
  y_lab <- function(lab) if (is_primary) lab else ""

  common_theme <- theme(text = element_text(size = 8), strip.background = element_blank())

  p <- list()
  p$height <- save_metric_plot(data = data_subset,
                               metric_col = "tree_height_prop_error",
                               y_label = y_lab("Tree Height\nProportional Error"),
                               show_x_labels = FALSE,
                               add_reference_line = TRUE,
                               reference_value = 0,
                               facet_col = facet_column,
                               ncol = 1) +
    common_theme + theme(strip.text = element_text(size = 8))

  p$length <- save_metric_plot(data = data_subset,
                               metric_col = "tree_length_prop_error",
                               y_label = y_lab("Tree Length\nProportional Error"),
                               show_x_labels = FALSE,
                               add_reference_line = TRUE,
                               reference_value = 0,
                               facet_col = facet_column,
                               ncol = 1) +
    common_theme + theme(strip.text = element_blank())

  p$rf <- save_metric_plot(data = data_subset,
                           metric_col = "rf_distance",
                           y_label = y_lab("RF Distance"),
                           show_x_labels = FALSE,
                           add_reference_line = TRUE,
                           reference_value = 0,
                           facet_col = facet_column,
                           ncol = 1) +
    common_theme + theme(strip.text = element_blank())

  p$mrca <- save_metric_plot(data = data_subset,
                             metric_col = "mrca_cell_type_accuracy",
                             y_label = y_lab("Ancestral Cell\nType Accuracy (%)"),
                             show_x_labels = TRUE,
                             add_reference_line = TRUE,
                             reference_value = 100,
                             facet_col = facet_column,
                             ncol = 1) +
    common_theme + theme(strip.text = element_blank())

  if (length(p) >= 4) {
    combined_plot <- p$height / p$length / p$rf / p$mrca +
      plot_layout(heights = c(1, 1, 1, 1)) &
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(1, 2, 1, 2)
      )
  }
  return(combined_plot)
}

# -------------------------- Execution --------------------------
if (plot_type == "main") {
  cat("Generating Main Figure (Selection)...\n")
  data_sel <- summary_data %>% filter(grepl("sel$", config))

  if (nrow(data_sel) > 0) {
    main_fig <- create_main_panels(data_sel %>% filter(grepl("tltt", simulation_name)), TRUE) |
      create_main_panels(data_sel %>% filter(grepl("gc_reentry", simulation_name)), FALSE)

    ggsave(file.path(pub_plots_dir, "main_figure_1to1_sel.pdf"), main_fig, width = 2.5, height = 5)
    cat("✓ Main figure saved\n")
  } else {
    stop("Error: No selection data found for main figure")
  }

} else if (plot_type == "supp") {
  cat("Generating Supplementary Figure (Neutral)...\n")
  data_neu <- summary_data %>% filter(grepl("neu$", config))

  if (nrow(data_neu) > 0) {
    supp_fig <- create_main_panels(data_neu %>% filter(grepl("tltt", simulation_name)), TRUE) |
      create_main_panels(data_neu %>% filter(grepl("gc_reentry", simulation_name)), FALSE)

    ggsave(file.path(supp_plots_dir, "supp_figure_1to1_neu.pdf"), supp_fig, width = 8, height = 7)
    cat("✓ Supplementary figure saved\n")
  } else {
    cat("Warning: No neutral data found. Skipping supplementary figure\n")
  }
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