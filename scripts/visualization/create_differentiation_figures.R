#!/usr/bin/env Rscript
# Create boxplots for cell type differentiation timing

suppressMessages({
  library(dowser)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(treeio)
})

source("scripts/utils/phylo_utilities.R")

# -------------------------- Command line Argument parsing --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript create_differentiation_figures.R job_list_csv array_task_id tree_analysis_dir pub_plots_dir\n")
  quit(status = 1)
}

job_list_file <- args[1]
array_task_id <- as.numeric(args[2])
tree_analysis_dir <- args[3]
pub_plots_dir <- args[4]

# Read job data
job_data <- read.csv(job_list_file, stringsAsFactors = FALSE)
job_row <- job_data[array_task_id, ]

config_name <- job_row$config_name
template_name <- job_row$template_name
true_tree_file <- job_row$true_tree_file
beast_tree_files <- strsplit(job_row$beast_tree_files, ";")[[1]]
analysis_type <- job_row$analysis_type

# Setup directories
plots_base_dir <- file.path(tree_analysis_dir, "plots")
dir.create(plots_base_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Differentiation Timing Job", array_task_id, "===\n")
cat("Config:", config_name, "\n")
cat("Template:", template_name, "\n")

# -------------------------- Helper functions --------------------------
# Helper function to normalize cell type names
normalize_tip <- function(x) {
  x <- as.character(x)
  tmp <- tolower(stringr::str_replace_all(x, "[ _]", ""))
  case_when(
    tmp %in% c("memorybcell", "memory_b_cell") ~ "MBC",
    tmp %in% c("plasmacell", "plasma_cell") ~ "PC",
    tmp %in% c("gcbcell", "gc", "default", "gc_b_cell") ~ "GC",
    TRUE ~ x
  )
}

# Find differentiation points by walking up the tree
getDiffPoint <- function(tree, node) {
  type <- filter(tree@data, node == !!node)$location
  edge <- tree@phylo$edge[tree@phylo$edge[, 2] == node, ]

  if (length(edge) == 0) {
    return(data.frame(diffnode = node, type = "root",
                      height = filter(tree@data, node == !!node)$height))
  }

  parent <- edge[1]
  parent_type <- filter(tree@data, node == parent)$location
  parent_height <- filter(tree@data, node == parent)$height

  if (parent_type == type) {
    getDiffPoint(tree, parent)
  } else {
    data.frame(diffnode = parent, type = parent_type, height = parent_height)
  }
}

# Get differentiation points for all tips
getDiffPoints <- function(tree) {
  diffpoints <- data.frame()
  for (tip in tree@phylo$tip.label) {
    tip_node <- which(tree@phylo$tip.label == tip)
    tip_data <- filter(tree@data, node == tip_node)
    diff_point <- getDiffPoint(tree, tip_node)

    temp <- data.frame(tip = tip, tip_type = tip_data$location, tip_height = tip_data$height)
    diffpoints <- rbind(diffpoints, cbind(temp, diff_point))
  }

  diffpoints$height <- as.numeric(diffpoints$height)
  diffpoints$tip_height <- as.numeric(diffpoints$tip_height)
  diffpoints
}

# Load true timing data from tree files
load_true_timing <- function(true_tree_file) {
  if (!file.exists(true_tree_file)) {
    cat("Warning: True tree file not found:", true_tree_file, "\n")
    return(data.frame())
  }

  cat("Loading true timing data from:", true_tree_file, "\n")

  tryCatch({
    true_trees <- treeio::read.beast(true_tree_file)
    true_data <- data.frame()

    for (clone_id in seq_along(true_trees)) {
      tree <- true_trees[[clone_id]]

      # Extract time_of_differentiation data
      if ("time_of_differentiation" %in% colnames(tree@data)) {
        clone_data <- tree@data %>%
          filter(!is.na(time_of_differentiation),
                 time_of_differentiation != "None",
                 time_of_differentiation != "") %>%
          mutate(clone_id = as.character(clone_id),
                 tip_type = normalize_tip(celltype),
                 tree_height = max(as.numeric(generation)) - 1,
                 # Convert time_of_differentiation to numeric
                 rel_differentiation_time = as.numeric(time_of_differentiation) / tree_height,
                 panel = "True") %>%
          filter(tip_type %in% c("MBC", "PC")) %>%
          select(clone_id, tip_type, rel_differentiation_time, panel)

        true_data <- rbind(true_data, clone_data)
      }
    }

    cat("Loaded true timing data for", length(unique(true_data$clone_id)), "clones\n")
    return(true_data)

  }, error = function(e) {
    cat("Error loading true timing data:", e$message, "\n")
    return(data.frame())
  })
}

# -------------------------- Main plotting function --------------------------
create_differentiation_timing_plots <- function() {

  cat("Processing BEAST tree files for template:", template_name, "\n")

  # Map clone IDs to beast files
  beast_tree_file_map <- list()
  clone_ids <- c()

  for (beast_tree_file in beast_tree_files) {
    extraction_result <- extract_template_and_clone(beast_tree_file)
    clone_id <- extraction_result$clone_id
    beast_tree_file_map[[clone_id]] <- beast_tree_file
    clone_ids <- c(clone_ids, clone_id)
  }

  clone_ids <- unique(clone_ids[order(as.numeric(clone_ids))])
  cat("Found", length(clone_ids), "clones:", paste(clone_ids, collapse = ", "), "\n")

  # Process each clone
  all_diffs <- list()
  tree_heights <- c()

  for (clone_id in clone_ids) {
    beast_tree_file <- beast_tree_file_map[[clone_id]]

    cat("Processing clone", clone_id, "...\n")

    tryCatch({
      beast_tree <- read.beast(beast_tree_file)
      diffs <- getDiffPoints(beast_tree)
      diffs$clone_id <- clone_id
      all_diffs[[clone_id]] <- diffs

      tree_height <- max(as.numeric(beast_tree@data$height), na.rm = TRUE)
      tree_heights[clone_id] <- tree_height

      cat("  Successfully processed clone", clone_id, "with tree height", tree_height, "\n")
    }, error = function(e) {
      cat("  Error processing clone", clone_id, ":", e$message, "\n")
    })
  }

  # Combine data
  estimated_diffs <- do.call(rbind, all_diffs)

  # Calculate relative heights and clean data
  estimated_diffs <- estimated_diffs %>%
    mutate(
      tree_height = tree_heights[as.character(clone_id)],
      relative_height = height / tree_height,
      tip_type = normalize_tip(tip_type),
      differentiation_time = - height,
      rel_differentiation_time = - relative_height + 1
    ) %>%
    filter(tip_type %in% c("MBC", "PC"))

  cat("Processed", length(unique(estimated_diffs$clone_id)), "clones successfully\n")

  # Load true data
  true_data <- load_true_timing(true_tree_file)

  # Create plots
  colors <- c("MBC" = "#0173B2", "PC" = "#E69F00")

  # Plot 1: Estimated relative timing
  p1 <- ggplot(estimated_diffs, aes(x = rel_differentiation_time, y = tip_type)) +
    geom_boxplot(aes(fill = tip_type), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = tip_type), width = 0, height = 0.1, size = 0.5) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      x = "Relative Differentiation Time (0 = root, 1 = present)",
      y = "Cell Type"
    ) +
    theme(legend.position = "none")

  # Plot 2: By clone
  p2 <- ggplot(estimated_diffs, aes(x = rel_differentiation_time, y = tip_type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = tip_type), width = 0, height = 0.1) +
    scale_color_manual(values = colors) +
    facet_wrap(~clone_id) +
    theme_bw() +
    labs(
      title = paste("Differentiation Timing by Clone -", template_name),
      subtitle = paste("Config:", config_name),
      x = "Differentiation Time (before present)",
      y = "Cell Type"
    ) +
    theme(legend.position = "none")

  plots_to_save <- list(p1, p2)

  # Plot 3: Estimated vs True
  if (nrow(true_data) > 0) {
    estimated_for_comparison <- estimated_diffs %>%
      select(clone_id, tip_type, rel_differentiation_time) %>%
      mutate(panel = "Estimated")

    combined_data <- rbind(
      estimated_for_comparison,
      true_data %>% select(clone_id, tip_type, rel_differentiation_time, panel)
    ) %>%
      mutate(
        panel = factor(panel, levels = c("True", "Estimated")),
        tip_type = factor(tip_type, levels = c("PC", "MBC"))
      )

    p3 <- ggplot(combined_data, aes(x = rel_differentiation_time, y = tip_type)) +
      geom_boxplot(aes(fill = tip_type), outlier.shape = NA) +
      geom_jitter(width = 0, height = 0.1, size = 0.025) +
      facet_grid(rows = vars(panel)) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      theme_bw() +
      labs(
        x = "Relative Differentiation Time (0 = root, 1 = present)",
        y = "Cell Type"
      ) +
      theme(
        legend.position = "none",
        strip.placement = "outside"
      )

    plots_to_save <- list(p1, p2, p3)

    # Save compact version
    compact_file <- file.path(pub_plots_dir, paste0("compact_differentiation_timing.pdf"))
    save_compact_plot(p3, compact_file, width = 2.5, height = 1.7)

  }

  # Save all plots using the utility function
  plot_file <- file.path(plots_base_dir, paste0(config_name, "_all_differentiation_timing.pdf"))
  save_plots_to_pdf(plots_to_save, plot_file, width = 12, height = 8)

  # Summary statistics
  summary_stats <- estimated_diffs %>%
    mutate(relative_time = 1 - relative_height) %>%
    group_by(tip_type) %>%
    summarize(
      n_observations = n(),
      n_clones = n_distinct(clone_id),
      mean_rel_time = mean(relative_time, na.rm = TRUE),
      median_rel_time = median(relative_time, na.rm = TRUE),
      sd_rel_time = sd(relative_time, na.rm = TRUE),
      .groups = 'drop'
    )

  summary_file <- file.path(plots_base_dir,
                            paste0(config_name, "_", template_name, "_timing_summary.csv"))
  write_csv(summary_stats, summary_file)

  cat("Plots saved to:", plot_file, "\n")
  cat("Summary saved to:", summary_file, "\n")

  return(list(estimated = estimated_diffs, true = true_data, summary = summary_stats))
}

# -------------------------- Main execution --------------------------
create_differentiation_timing_plots()
cat("Differentiation timing job", array_task_id, "completed\n")