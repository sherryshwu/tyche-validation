#!/usr/bin/env Rscript
# Load required packages
if (!exists("packages_loaded")) {
  suppressMessages({
    library(treeio)
    library(ggtree)
    library(ape)
    library(tidyverse)
    library(tidyr)
    library(phangorn)
    library(patchwork)
    library(readr)
  })
  packages_loaded <- TRUE
}

# -------------------------- Utility Functions --------------------------
# Build vector of parameters to ignore during convergence checking
build_ignore_vector <- function(template_id, model_type, analysis_scope) {
  is_fixed   <- grepl("FixedClockRates", template_id)
  is_relaxed <- grepl("UCRelaxedClock", template_id)
  is_3state  <- grepl("3state", template_id)

  ignore <- c("traitfrequencies", "freqParameter", "rateIndicator")

  # Analysis-specific ignore parameters
  if (analysis_scope == "differentiation_analysis" && is_3state) {
    ignore <- c(ignore, "traitRates", "clockRate", "traitClockRate")
  }

  if (is_fixed) {
    ignore <- c(ignore, "typeLinkedRates")
    if (identical(model_type, "competing_models")) {
      ignore <- c(ignore, "^geneticClockRate")
    }
  }
  if (is_relaxed) {
    ignore <- c(ignore, "rateCategories")
  }

  unique(ignore)
}

# Ensure neutral evolution fields are present in data
ensure_neutral_fields <- function(df) {
  if (!"productive" %in% names(df))                df$productive <- TRUE
  if (!"v_call" %in% names(df))                    df$v_call <- "IGHV1-1*01"
  if (!"d_call" %in% names(df))                    df$d_call <- "IGHD1-1*01"
  if (!"j_call" %in% names(df))                    df$j_call <- "IGHJ1*01"
  if (!"rev_comp" %in% names(df))                  df$rev_comp <- FALSE
  if (!"germline_alignment_d_mask" %in% names(df)) df$germline_alignment_d_mask <- df$germline_alignment
  df
}

# Standardized model naming abbreviation
get_model_short_name <- function(model_name) {
  case_when(
    grepl("ExpectedOccupancy.*Est", model_name) ~ "EO_Est",
    grepl("ExpectedOccupancy.*Fixed", model_name) ~ "EO_Fixed",
    grepl("InstantSwitch.*Est", model_name) ~ "IS_Est",
    grepl("MixedSwitch.*Est", model_name) ~ "MS_Est",
    grepl("MixedSwitch.*Fixed", model_name) ~ "MS_Fixed",
    grepl("StrictClock.*AncestralReconstruction", model_name) ~ "SC_AR",
    grepl("UCRelaxedClock.*AncestralReconstruction", model_name) ~ "UCLD_AR",
    TRUE ~ model_name
  )
}

# Extract template name and clone ID from filename
extract_template_and_clone <- function(filename) {
  base_name <- basename(tools::file_path_sans_ext(filename))
  pattern <- "^(.+)_([0-9]+)_tree_with_trait.*$"

  if (grepl(pattern, base_name)) {
    matches <- regmatches(base_name, regexec(pattern, base_name))[[1]]
    return(list(template_name = matches[2], clone_id = matches[3]))
  }
  return(list(template_name = "Unknown", clone_id = "1"))
}

# Helper function to read ESS from BEAST log file
read_ess_from_log <- function(beast_tree_file) {
  # Construct log file path
  log_tsv_file <- gsub("_tree_with_trait\\.tree$", "_log.tsv", beast_tree_file)

  if (!file.exists(log_tsv_file)) {
    cat("  Warning: Log file not found:", basename(log_tsv_file), "\n")
    return(NA)
  }
  # Read TSV file
  lines <- readLines(log_tsv_file, warn = FALSE)
  tree_height_line <- grep("^TreeHeight\\s+", lines, value = TRUE)
  tree_height_ess <- as.numeric(unlist(strsplit(tree_height_line[1], "\\s+"))[9])

  return(tree_height_ess)
}

# Create standardized directory structure
create_plot_directories <- function(base_dir, plot_types) {
  dirs <- list()
  for (plot_type in plot_types) {
    dir_path <- file.path(base_dir, plot_type)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    dirs[[plot_type]] <- dir_path
  }
  return(dirs)
}

# -------------------------- Plotting Utility Functions --------------------------
# Standard color schemes
get_location_colors <- function(analysis_type = "main_analysis") {
  if (analysis_type %in% c("main_analysis", "sub_analysis")) {
    types <- RColorBrewer::brewer.pal(3, "Set1")
    names(types) <- c("germinal_center", "other", "germinal_center+other")
  } else if (analysis_type == "differentiation_analysis") {
    types <- RColorBrewer::brewer.pal(3, "Set2")
    names(types) <- c("gc_b_cell", "plasma_cell", "memory_b_cell")
  } else {
    stop("Unknown analysis_type: use main_analysis, sub_analysis or differentiation_analysis")
  }
  return(types)
}

# Always uses 2-value scheme for non-beast trees
get_standard_tree_colors <- function() {
  types <- RColorBrewer::brewer.pal(3, "Set1")
  names(types) <- c("germinal_center", "other", "germinal_center+other")
  return(types)
}

# Standardize location labels between different datasets
standardize_location_labels <- function(tree_data) {
  if ("location" %in% colnames(tree_data)) {
    tree_data$location <- case_when(
      tree_data$location == "GC" ~ "germinal_center",
      tree_data$location == "GC+other" ~ "germinal_center+other",
      TRUE ~ tree_data$location
    )
  }
  return(tree_data)
}

save_metric_plot <- function(metric_col, y_label, title = NULL,
                             show_x_labels = TRUE, add_reference_line = FALSE,
                             reference_value = NULL, facet_col = "facet_label_with_sim",
                             nrow = 1, ncol = NULL, facet_scales = "fixed") {
  plot_data <- summary_data
  plot_data$metric_value <- plot_data[[metric_col]]

  p <- ggplot(plot_data, aes(x = model_display, y = metric_value, fill = model_display)) +
    geom_boxplot(alpha = 1, outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.2, height = 0, size = 0.2) +
    facet_wrap(as.formula(paste("~", facet_col)), scales = facet_scales, nrow = nrow, ncol = ncol) +
    scale_fill_brewer(palette = "Pastel1") +
    labs(x = "", y = y_label, fill = "Model", title = title) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_line(linewidth = 0.2),
          panel.grid.minor = element_line(linewidth = 0.1),
          axis.line = element_blank(),
          panel.border = element_rect(
            color = "black",
            linewidth = 0.3,
            fill = NA
          ))

  # Add reference line if specified
  if (add_reference_line && !is.null(reference_value)) {
    p <- p + geom_hline(yintercept = reference_value, linetype = "dashed", linewidth = 0.3)
  }

  # Handle x-axis labels
  if (!show_x_labels) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(p)
}

save_compact_plot <- function(plot, filename = NULL, width = 3, height = 1.5,
                              legend.position = "none", box_size = 0.15,
                              line_size = 0.5, text_size = 8, font_family = "Helvetica",
                              dpi = 300, extra_theme = NULL,
                              remove_x_labels = FALSE) {
  base_theme <- theme(
    text = element_text(family = font_family, size = text_size, color = "black"),
    plot.title = element_text(size = text_size, hjust = 0.5),
    axis.title = element_text(size = text_size),
    axis.text.x = element_text(size = text_size, color = "black", hjust = 0.5),
    axis.text.y = element_text(size = text_size, color = "black"),
    axis.ticks = element_line(linewidth = line_size, color = "black"),
    panel.border = element_rect(color = "black", linewidth = line_size, fill = NA),
    strip.text = element_text(size = 8, face = "plain", margin = margin(2, 0, 4, 0)),
    strip.background = element_blank(),
    legend.position = legend.position,
    legend.key.size = unit(0.35, "lines"),
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = text_size),
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

  if (remove_x_labels) {
    base_theme <- base_theme + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  }

  plot <- plot + base_theme

  if (!is.null(extra_theme)) {
    plot <- plot + extra_theme
  }

  plot$layers <- lapply(plot$layers, function(layer) {
    if (inherits(layer$geom, "GeomBoxplot")) {
      layer$aes_params$linewidth <- box_size
    }
    return(layer)
  })

  ggsave(filename, plot, width = width, height = height, units = "in", dpi = dpi)
  return(plot)
}

# Basic plot saving to PDF (for non-compact plots)
save_plots_to_pdf <- function(plot_list, filename, width = 12, height = 8) {
  if (length(plot_list) == 0) {
    warning("No plots to save")
    return(FALSE)
  }

  pdf(filename, width = width, height = height)
  for (plot in plot_list) {
    if (inherits(plot, "ggplot")) {
      print(plot)
    }
  }
  dev.off()

  cat("Saved", length(plot_list), "plots to", basename(filename), "\n")
  TRUE
}