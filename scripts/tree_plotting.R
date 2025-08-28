#!/usr/bin/env Rscript
source("scripts/tree_functions.R")
suppressMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(dplyr)
  library(stringr)
})
# -------------------------- Command line Argument parsing --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript tree_plotting_runner.R job_list_csv array_task_id tree_analysis_dir\n")
  quit(status = 1)
}

job_list_file <- args[1]
array_task_id <- as.numeric(args[2])
tree_analysis_dir <- args[3]

# Read job data
job_data <- read.csv(job_list_file, stringsAsFactors = FALSE)
job_row <- job_data[array_task_id, ]

config_name <- job_row$config_name
template_name <- job_row$template_name
true_tree_file <- job_row$true_tree_file
beast_files <- strsplit(job_row$beast_files, ";")[[1]]
model_type <- job_row$model_type

# Setup directories
plots_base_dir <- file.path(tree_analysis_dir, "plots", "tree_plots")
dir.create(plots_base_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Tree Plotting Job", array_task_id, "===\n")
cat("Config:", config_name, "\n")
cat("Template:", template_name, "\n")
cat("Model type:", model_type, "\n")

# -------------------------- Plotting functions --------------------------
create_all_tree_plots <- function(beast_tree_files, true_tree_file, plots_dir,
                                  config_name, template_name) {

  cat("=== Creating comprehensive tree plots ===\n")

  # Create plot directories
  plot_dirs <- create_plot_directories(plots_base_dir, c(
    "time_trees", "cophylo_plots", "genetic_distance_trees", "combined_plots"
  ))

  location_colors <- get_location_colors()

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
  beast_file_map <- list()
  clone_ids <- c()

  for (beast_tree_file in beast_tree_files) {
    extraction_result <- extract_template_and_clone(beast_tree_file)
    clone_id <- extraction_result$clone_id
    beast_file_map[[clone_id]] <- beast_tree_file
    clone_ids <- c(clone_ids, clone_id)
  }

  clone_ids <- unique(clone_ids[order(as.numeric(clone_ids))])

  # Initialize plot collections
  all_plots <- list(
    comparison = list(),
    cophylo = list(),
    genetic = list(),
    combined = list()
  )

  # Generate plots for each clone
  for (clone_id in clone_ids) {
    tryCatch({
      cat("Processing clone", clone_id, "...\n")

      beast_tree_file <- beast_file_map[[clone_id]]
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
      beast_plot <- plot_time_tree(beast_tree, paste("BEAST - Clone", clone_id), location_colors)
      true_time_plot <- plot_time_tree(true_time_tree, paste("True Time - Clone", clone_id), location_colors)

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
                                              template_name, config_name, location_colors)
      if (!is.null(cophylo_plot)) {
        all_plots$cophylo[[clone_id]] <- cophylo_plot
      }

      # 3. Genetic distance tree
      true_genetic_tree <- true_genetic_trees[[as.numeric(clone_id)]]
      if (!is.null(true_genetic_tree)) {
        genetic_plot <- plot_genetic_tree(true_genetic_tree, clone_id, config_name, location_colors)
        if (!is.null(genetic_plot)) {
          all_plots$genetic[[clone_id]] <- genetic_plot
        }

        # 4. Combined plot
        genetic_plot_small <- plot_genetic_tree(true_genetic_tree, clone_id,
                                                config_name, location_colors, compact = TRUE)

        combined_plot <- beast_plot + true_time_plot + genetic_plot_small +
          plot_layout(ncol = 3, guides = "collect") &
          theme(legend.position = "bottom")

        combined_plot <- combined_plot +
          plot_annotation(
            title = paste("All trees - Clone", clone_id),
            subtitle = paste("BEAST | True Time | True Genetic | Template:", template_name)
          )

        all_plots$combined[[clone_id]] <- combined_plot
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

plot_genetic_tree <- function(tree, clone_id, config_name, location_colors, compact = FALSE) {
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

create_cophylogeny_plot <- function(beast_tree, true_tree, clone_id, template_name, config_name, location_colors) {
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
      scale_color_manual(values = location_colors) +
      theme_tree()

    p2 <- ggtree(true_tree_rotated, ladderize = FALSE) +
      geom_tippoint(aes(color = location), size = 1.5) +
      geom_nodepoint(aes(color = location), size = 1) +
      scale_color_manual(values = location_colors) +
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
      scale_color_manual(values = location_colors)

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
         name = "combined_trees", width = 18, height = 8)
  )

  for (config in plot_configs) {
    if (length(config$plots) > 0) {
      filename <- paste0(config_name, "_", config$name, "_", template_name, "_all_clones.pdf")
      save_plots_to_pdf(config$plots, file.path(config$dir, filename), config$width, config$height)
      cat("✓ Saved", config$name, "\n")
    }
  }
}

# -------------------------- Main execution --------------------------
create_all_tree_plots(
  beast_tree_files = beast_files,
  true_tree_file = true_tree_file,
  plots_dir = plots_base_dir,
  config_name = config_name,
  template_name = template_name
)

cat("✓ Tree plotting job", array_task_id, "completed\n")