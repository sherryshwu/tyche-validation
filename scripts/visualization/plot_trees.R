#!/usr/bin/env Rscript
source("scripts/utils/phylo_utilities.R")
suppressMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(dplyr)
  library(stringr)
  library(cowplot)
  library(RColorBrewer)
  library(ape)
})
# -------------------------- Command line Argument parsing --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript plot_trees.R job_list_csv array_task_id tree_analysis_dir analysis_type\n")
  quit(status = 1)
}

job_list_file <- args[1]
array_task_id <- as.numeric(args[2])
tree_analysis_dir <- args[3]
analysis_type <- args[4]

# Read job data
job_data <- read.csv(job_list_file, stringsAsFactors = FALSE)
job_row <- job_data[array_task_id, ]

analysis_type <- job_row$analysis_type
config_name <- job_row$config_name
template_name <- job_row$template_name
true_tree_file <- job_row$true_tree_file
beast_tree_files <- strsplit(job_row$beast_tree_files, ";")[[1]]

# Setup directories
plots_base_dir <- file.path(tree_analysis_dir, "plots", "tree_plots")
dir.create(plots_base_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Tree Plotting Job", array_task_id, "===\n")
cat("Analysis type:", analysis_type, "\n")
cat("Config:", config_name, "\n")
cat("Template:", template_name, "\n")

# -------------------------- Plotting functions --------------------------
create_all_tree_plots <- function(beast_tree_files, true_tree_file, plots_dir,
                                  config_name, template_name, analysis_type) {

  cat("=== Creating comprehensive tree plots ===\n")

  # Create plot directories
  plot_dirs <- create_plot_directories(plots_base_dir, c(
    "time_trees", "cophylo_plots", "genetic_distance_trees", "combined_plots", "four_tree_plots"
  ))

  # BEAST trees use analysis-specific colors
  beast_location_colors <- get_location_colors(analysis_type)

  # True trees and genetic trees always use standard 2-value colors
  standard_location_colors <- get_standard_tree_colors()

  # Load trees
  cat("Loading trees...\n")
  tryCatch({
    # Load true time trees
    true_time_trees <- read.beast(true_tree_file)
    cat("✓ Loaded", length(true_time_trees), "time trees\n")

    # Load true genetic distance trees
    genetic_tree_file <- gsub("all_simplified_time_trees\\.nex$",
                              "all_simplified_trees.nex", true_tree_file)

    true_genetic_trees <- NULL
    if (file.exists(genetic_tree_file)) {
      true_genetic_trees <- read.beast(genetic_tree_file)
      cat("✓ Loaded genetic distance trees\n")
    }
  }, error = function(e) {
    stop("Failed to load tree files: ", e$message)
  })

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

  # Initialize plot collections
  all_plots <- list(
    comparison = list(),
    cophylo = list(),
    genetic = list(),
    combined = list(),
    four_trees = list()
  )

  # Generate plots for each clone
  for (clone_id in clone_ids) {
    tryCatch({
      cat("Processing clone", clone_id, "...\n")

      beast_tree_file <- beast_tree_file_map[[clone_id]]
      if (is.null(beast_tree_file)) {
        cat(" No BEAST tree file for clone", clone_id, "\n")
        next
      }

      beast_tree <- read.beast(beast_tree_file)
      true_time_tree <- true_time_trees[[as.numeric(clone_id)]]

      if (is.null(true_time_tree)) {
        cat(" True time tree not found for clone", clone_id, "\n")
        next
      }

      # 1. Time tree comparison
      beast_plot <- plot_time_tree(beast_tree, paste("BEAST - Clone", clone_id), beast_location_colors)
      true_time_plot <- plot_time_tree(true_time_tree, paste("True Time - Clone", clone_id), standard_location_colors)

      comparison_plot <- beast_plot + true_time_plot +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

      comparison_plot <- comparison_plot +
        plot_annotation(
          title = paste("Time Tree Comparison - Clone", clone_id),
          subtitle = paste("Template:", template_name, "| Config:", config_name)
        )

      all_plots$comparison[[clone_id]] <- comparison_plot

      # 2. Cophylogeny plot
      cophylo_plot <- create_cophylogeny_plot(beast_tree, true_time_tree, clone_id,
                                              template_name, config_name,
                                              beast_location_colors, standard_location_colors)
      if (!is.null(cophylo_plot)) {
        all_plots$cophylo[[clone_id]] <- cophylo_plot
      }

      # 3. Genetic distance tree
      true_genetic_tree <- true_genetic_trees[[as.numeric(clone_id)]]
      if (!is.null(true_genetic_tree)) {
        genetic_plot <- plot_genetic_tree(true_genetic_tree, clone_id, config_name)
        if (!is.null(genetic_plot)) {
          all_plots$genetic[[clone_id]] <- genetic_plot
        }

        # 4. Combined plot
        genetic_plot_small <- plot_genetic_tree(true_genetic_tree, clone_id, config_name, compact = TRUE)

        combined_plot <- beast_plot + true_time_plot + genetic_plot_small +
          plot_layout(ncol = 3, guides = "collect") &
          theme(legend.position = "bottom")

        combined_plot <- combined_plot +
          plot_annotation(
            title = paste("All trees - Clone", clone_id),
            subtitle = paste("BEAST | True Time | True Genetic | Template:", template_name)
          )

        all_plots$combined[[clone_id]] <- combined_plot

        # 5. Four-tree scaled plot
        four_tree_plot <- create_four_tree_plot(beast_tree_file, true_time_tree, true_genetic_tree, 
                                                clone_id, template_name, config_name, 
                                                beast_location_colors, standard_location_colors)
        if (!is.null(four_tree_plot)) {
          all_plots$four_trees[[clone_id]] <- four_tree_plot
        }
      }
      cat("✓ Plots created for clone", clone_id, "\n")

    }, error = function(e) {
      cat("✗ ERROR for clone", clone_id, ":", e$message, "\n")
    })
  }

  # Save all plots
  save_all_plots(all_plots, plot_dirs, template_name, config_name)

  cat("All tree plotting completed\n")
}

plot_time_tree <- function(tree, title, location_colors) {
  tryCatch({
    # Standardize location labels
    tree@data <- standardize_location_labels(tree@data)

    p <- ggtree(tree) +
      geom_point(aes(color = if ("location" %in% colnames(tree@data)) location else NULL), size = 2) +
      geom_treescale() +
      theme_tree() +
      scale_color_manual(values = location_colors) +
      labs(title = title, color = "Location") +
      theme(plot.title = element_text(size = 11, hjust = 0.5),
            legend.position = "bottom")

    return(p)

  }, error = function(e) {
    cat("Error occurred when plotting time trees!")
  })
}

plot_genetic_tree <- function(tree, clone_id, config_name, compact = FALSE) {
  location_colors <- get_standard_tree_colors()
  tryCatch({
    if (!"location" %in% colnames(tree@data)) {
      return(NULL)
    }

    p <- ggtree(tree) +
      geom_tippoint(aes(color = location), size = if (compact) 1 else 2) +
      scale_color_manual(values = location_colors) +
      theme_tree() +
      labs(
        title = if (compact) paste("Genetic - Clone", clone_id) else "Genetic Distance Tree",
        subtitle = if (!compact) paste("Config:", config_name, "- Clone", clone_id) else NULL,
        color = "Location"
      ) +
      theme(
        plot.title = element_text(size = if (compact) 10 else 11, hjust = 0.5),
        legend.position = if (compact) "none" else "bottom"
      )

    return(p)

  }, error = function(e) {
    cat("Error occurred when plotting genetic distance trees!")
  })
}

create_cophylogeny_plot <- function(beast_tree, true_tree, clone_id, template_name, config_name, beast_location_colors, true_location_colors) {
  tryCatch({
    # Convert to phylo objects
    beast_phylo <- as.phylo(beast_tree)
    true_phylo <- as.phylo(true_tree)

    # Standardize tip labels
    beast_phylo$tip.label <- gsub("_heavy$", "", beast_phylo$tip.label)
    common_tips <- intersect(beast_phylo$tip.label, true_phylo$tip.label)

    if (length(common_tips) == 0) {
      cat("  Warning: No common tip labels found for cophylogeny plot\n")
      return(NULL)
    }

    # Clean and prepare trees
    beast_tree_clean <- beast_tree
    beast_tree_clean@phylo$tip.label <- gsub("_heavy$", "", beast_tree_clean@phylo$tip.label)

    beast_plot_temp <- ggtree(beast_tree_clean)
    beast_order <- get_taxa_name(beast_plot_temp)
    beast_order <- beast_order[!is.na(beast_order) & beast_order != ""]

    tips_to_drop <- setdiff(true_tree@phylo$tip.label, beast_order)
    true_phylo_subset <- if (length(tips_to_drop) > 0) {
      drop.tip(true_tree@phylo, tips_to_drop)
    } else {
      true_tree@phylo
    }

    true_tree_subset <- true_tree
    true_tree_subset@phylo <- true_phylo_subset

    if (!ape::is.binary(true_tree_subset@phylo)) {
      true_tree_subset_binary <- true_tree_subset
      true_tree_subset_binary@phylo <- ape::multi2di(true_tree_subset@phylo, random = FALSE)
    } else {
      true_tree_subset_binary <- true_tree_subset
    }

    # Rotate trees for alignment
    true_tree_rotated <- true_tree_subset_binary
    true_tree_rotated@phylo <- ape::rotateConstr(true_tree_subset_binary@phylo, beast_order)
    beast_tree_rotated <- beast_tree_clean
    beast_tree_rotated@phylo <- ape::rotateConstr(beast_tree_clean@phylo, beast_order)

    # Create plots
    p1 <- ggtree(beast_tree_rotated, ladderize = FALSE) +
      geom_tippoint(aes(color = location), size = 1.5) +
      geom_nodepoint(aes(color = location), size = 1) +
      scale_color_manual(values = beast_location_colors) +
      theme_tree()

    p2 <- ggtree(true_tree_rotated, ladderize = FALSE) +
      geom_tippoint(aes(color = location), size = 1.5) +
      geom_nodepoint(aes(color = location), size = 1) +
      scale_color_manual(values = true_location_colors) +
      theme_tree()

    # Combine plots
    d1 <- p1$data
    d2 <- p2$data
    d2$x <- max(d2$x) - d2$x + max(d1$x) + 5

    pp <- p1 +
      geom_tree(data = d2) +
      geom_tippoint(data = d2, aes(color = location), size = 1.5) +
      geom_nodepoint(data = d2, aes(color = location), size = 1) +
      ggnewscale::new_scale_color() +
      scale_color_manual(values = true_location_colors)

    # Add connecting lines if possible
    if ("label" %in% colnames(d1) && "label" %in% colnames(d2)) {
      dd <- bind_rows(d1, d2) %>%
        filter(!is.na(label), label != "", label != "NA") %>%
        mutate(clean_label = gsub("_heavy$", "", label)) %>%
        filter(clean_label %in% common_tips) %>%
        filter(!is.na(clean_label), clean_label != "")

      if (nrow(dd) > 0) {
        pp <- pp +
          geom_line(aes(x, y, group = clean_label), data = dd,
                    color = 'lightblue', alpha = 0.6, linewidth = 0.3)
      }
    }

    final_plot <- pp +
      labs(title = paste("Cophylogeny Plot - Clone", clone_id),
           subtitle = paste("Template:", template_name, "| Config:", config_name),
           color = "Location") +
      theme(plot.title = element_text(size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            legend.position = "bottom")

    return(final_plot)

  }, error = function(e) {
    cat("  ✗ Error creating cophylogeny plot:", e$message, "\n")
    return(NULL)
  })
}

save_all_plots <- function(all_plots, plot_dirs, template_name, config_name) {
  plot_configs <- list(
    list(plots = all_plots$comparison, dir = plot_dirs$time_trees,
         name = "time_tree_comparison", width = 14, height = 8),
    list(plots = all_plots$cophylo, dir = plot_dirs$cophylo_plots,
         name = "cophylo", width = 16, height = 8),
    list(plots = all_plots$genetic, dir = plot_dirs$genetic_distance_trees,
         name = "genetic_distance_trees", width = 10, height = 8),
    list(plots = all_plots$combined, dir = plot_dirs$combined_plots,
         name = "combined_trees", width = 18, height = 8),
    list(plots = all_plots$four_trees, dir = plot_dirs$four_tree_plots,
         name = "four_tree_comparison", width = 8, height = 3.25)
  )

  for (config in plot_configs) {
    if (length(config$plots) > 0) {
      filename <- paste0(config_name, "_", config$name, "_", template_name, "_all_clones.pdf")
      save_plots_to_pdf(config$plots, file.path(config$dir, filename), config$width, config$height)
      cat("✓ Saved", config$name, "\n")
    }
  }
}

create_four_tree_plot <- function(beast_tree_file, true_time_tree, true_genetic_tree, 
                                  clone_id, template_name, config_name, 
                                  beast_location_colors, standard_location_colors) {
  tryCatch({
    # Find EO_Fixed and SC_AR models for this clone
    eo_fixed_file <- NULL
    sc_ar_file <- NULL

    # Extract the base path structure from the current beast_tree_file
    base_analysis_path <- sub("(.*rev).*", "\\1", beast_tree_file)
    config_dir <- basename(dirname(beast_tree_file))

    # Look for EO_Fixed (type-linked) tree
    eo_fixed_file <- file.path(base_analysis_path, "tyche_models", "beast_raw_output",
                               config_dir, paste0("ExpectedOccupancy_FixedClockRates_", clone_id, "_tree_with_trait.tree"))

    # Look for SC_AR (strict clock) tree
    sc_ar_file <- file.path(base_analysis_path, "competing_models", "beast_raw_output",
                            config_dir, paste0("StrictClock_AncestralReconstruction_", clone_id, "_tree_with_trait.tree"))

    cat("    Looking for EO_Fixed:", eo_fixed_file, "\n")
    cat("    Looking for SC_AR:", sc_ar_file, "\n")

    # Skip if we don't have both required models
    if (is.null(eo_fixed_file) || is.null(sc_ar_file)) {
      cat("    Missing EO_Fixed or SC_AR for clone", clone_id, "\n")
      return(NULL)
    }

    # Load BEAST trees
    eo_fixed_tree <- read.beast(eo_fixed_file)
    sc_ar_tree <- read.beast(sc_ar_file)

    # Calculate maximum tree height for scaling
    tree_height <- max(
      max(as.numeric(node.depth.edgelength(true_time_tree@phylo)) - 1),
      max(as.numeric(node.depth.edgelength(eo_fixed_tree@phylo))),
      max(as.numeric(node.depth.edgelength(sc_ar_tree@phylo)))
    )

    # Plot parameters
    tree_size <- 0.25
    node_size <- 1.5
    tip_size <- 1.5
    stroke_size <- 0.25

    # Create individual plots
    gdp <- create_scaled_tree_plot(true_genetic_tree, "GD Parsimony", standard_location_colors,
                                   tree_size, node_size, tip_size, stroke_size)

    timep <- create_scaled_tree_plot(true_time_tree, "True", standard_location_colors,
                                     tree_size, node_size, tip_size, stroke_size, 
                                     tree_height = tree_height)

    tlp <- create_scaled_tree_plot(eo_fixed_tree, "Type-linked", beast_location_colors,
                                   tree_size, node_size, tip_size, stroke_size, 
                                   tree_height = tree_height)

    slp <- create_scaled_tree_plot(sc_ar_tree, "Strict", beast_location_colors,
                                   tree_size, node_size, tip_size, stroke_size, 
                                   tree_height = tree_height)

    # Create legend
    legend_plot <- ggplot() +
      geom_point(data = data.frame(x = seq_along(standard_location_colors), y = 1,
                                   location = names(standard_location_colors)),
                 aes(x = x, y = y, fill = location), size = 3, shape = 21) +
      scale_fill_manual(values = standard_location_colors, name = "Location") +
      theme_void() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8))

    legend <- cowplot::get_legend(legend_plot)

    # Remove axes theme
    no_y_axis <- theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),
                       panel.border = element_blank(),
                       plot.title = element_text(hjust = 0.5))

    # Apply compact styling
    gdp_compact <- gdp + no_y_axis + theme(legend.position = "none")
    timep_compact <- timep + no_y_axis + theme(legend.position = "none")
    tlp_compact <- tlp + no_y_axis + theme(legend.position = "none")
    slp_compact <- slp + no_y_axis + theme(legend.position = "none")

    # Combine plots
    combined_plot <- plot_grid(gdp_compact, timep_compact, tlp_compact, slp_compact,
                               nrow = 1, align = "h")

    # Add legend
    final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.15))

    # Add title
    final_plot <- final_plot +
      plot_annotation(
        title = paste(template_name, "- Clone", clone_id)
      )

    return(final_plot)

  }, error = function(e) {
    cat("    ✗ Error creating four-tree plot:", e$message, "\n")
    return(NULL)
  })
}

create_scaled_tree_plot <- function(tree, title, location_colors, tree_size, node_size, tip_size, stroke_size, tree_height = NULL) {
  tryCatch({
    # Standardize location labels
    tree@data <- standardize_location_labels(tree@data)

    p <- ggtree(tree, size = tree_size) + 
      geom_nodepoint(aes(fill = location), size = node_size, pch = 21, stroke = stroke_size) +
      geom_tippoint(aes(fill = location), size = tip_size, pch = 21, stroke = stroke_size) +
      scale_fill_manual(values = location_colors) +
      labs(fill = "Location", title = title) +
      coord_cartesian(ylim = c(0, 100))

    # Add tree scales
    if (is.null(tree_height)) {
      # For genetic distance trees, use appropriate scale
      p <- p + geom_treescale(width = 15, linesize = 0.15, fontsize = 2, offset = 0.5)
    } else {
      # For time trees, use time-appropriate scale
      p <- p + geom_treescale(width = 100, linesize = 0.15, fontsize = 2, offset = 0.5)
    }

    # Apply xlim for time trees to ensure same scale
    if (!is.null(tree_height)) {
      p <- p + xlim(0, tree_height)
    }

    return(p)

  }, error = function(e) {
    cat("Error plotting tree:", e$message, "\n")
    return(ggplot() + theme_void() + labs(title = paste("Error:", title)))
  })
}

# -------------------------- Main execution --------------------------
create_all_tree_plots(
  beast_tree_files = beast_tree_files,
  true_tree_file = true_tree_file,
  plots_dir = plots_base_dir,
  config_name = config_name,
  template_name = template_name,
  analysis_type = analysis_type
)

cat("✓ Tree plotting job", array_task_id, "completed\n")