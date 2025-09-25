#!/usr/bin/env Rscript
# Load libraries
suppressMessages({
  library(treeio)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(dowser)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat("Usage: Rscript tree_analysis.R job_list_csv array_task_id output_dir simulation_name analysis_type rev_suffix\n")
  quit(status = 1)
}

job_list_file <- args[1]
array_task_id <- as.numeric(args[2])
output_dir <- args[3]
simulation_name <- args[4]  # tltt_08_20, gc_reentry_hunter, etc.
analysis_type <- args[5]    # main_analysis, differentiation_analysis, etc.
rev_suffix <- args[6]       # irrev, rev

# Read CSV and get the specific job
job_data <- read.csv(job_list_file)
job_row <- job_data[array_task_id, ]

config_name <- job_row$config_name
template_name <- job_row$template_name
true_tree_file <- job_row$true_tree_file
beast_tree_files <- strsplit(as.character(job_row$beast_tree_files), ";")[[1]]

cat("=== Tree Analysis Started ===\n")
cat("Simulation:", simulation_name, "\n")
cat("Analysis:", analysis_type, "\n") 
cat("Rev suffix:", rev_suffix, "\n")
cat("Config:", config_name, "\n")
cat("Template:", template_name, "\n")
cat("BEAST files:", length(beast_tree_files), "\n")

# Create output directory
config_output_dir <- file.path(output_dir, "results", config_name)
dir.create(config_output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to read TreeHeight from BEAST log file
read_tree_height_from_log <- function(beast_tree_file) {
  # Construct log file path
  log_tsv_file <- gsub("_tree_with_trait\\.tree$", "_log.tsv", beast_tree_file)

  if (!file.exists(log_tsv_file)) {
    cat("  Warning: Log file not found:", basename(log_tsv_file), "\n")
    return(NA)
  }
  # Read TSV file
  lines <- readLines(log_tsv_file, warn = FALSE)
  tree_height_line <- grep("^TreeHeight\\s+", lines, value = TRUE)
  tree_height_mean <- unlist(strsplit(tree_height_line[1], "\\s+"))[2]

  return(tree_height_mean)
}

# Helper function for MRCA analysis
standardize_mrca_locations <- function(mrca_data) {
  if ("mrca_location" %in% colnames(mrca_data)) {
    mrca_data$mrca_location <- dplyr::case_when(
      mrca_data$mrca_location == "GC" ~ "germinal_center",
      mrca_data$mrca_location == "GC+other" ~ "germinal_center+other",
      TRUE ~ mrca_data$mrca_location
    )
  }
  mrca_data
}

extract_mrca_locations <- function(tree, clean_tip_label = FALSE) {
  tree_phylo <- as.phylo(tree)
  tips <- tree_phylo$tip.label

  if (clean_tip_label) {
    tips <- gsub("_heavy$", "", tips)
    tree_phylo$tip.label <- tips
  }

  n_tips <- length(tips)

  if ("location" %in% colnames(tree@data)) {
    node_data <- tree@data %>%
      select(location, node) %>%
      mutate(node_num = as.numeric(node))

    missing_locations <- sum(is.na(node_data$location))
    if (missing_locations > 0) {
      cat("  Missing locations:", missing_locations, "\n")
    }
  } else {
    cat("  WARNING: No 'location' column found in tree data\n")
    return(data.frame())
  }

  mrca_data <- data.frame(
    tip1 = character(),
    tip2 = character(),
    mrca_node = numeric(),
    mrca_location = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:(n_tips - 1)) {
    for (j in (i + 1):n_tips) {
      mrca_node <- getMRCA(tree_phylo, c(tips[i], tips[j]))

      if (!is.null(mrca_node) && !is.na(mrca_node)) {
        node_idx <- which(node_data$node_num == mrca_node)

        if (length(node_idx) > 0) {
          mrca_location <- node_data$location[node_idx]
        } else {
          mrca_location <- NA
        }

        mrca_data <- rbind(mrca_data, data.frame(
          tip1 = tips[i],
          tip2 = tips[j],
          mrca_node = mrca_node,
          mrca_location = mrca_location,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(mrca_data)
}

# Load true trees
cat("Loading true trees...\n")
true_trees <- tryCatch({
  treeio::read.beast(true_tree_file)
}, error = function(e) {
  cat("ERROR: Failed to load true trees:", e$message, "\n")
  quit(status = 1)
})

cat("Loaded", length(true_trees), "true trees\n")

# Process each BEAST tree file
results_list <- list()
summary_data <- data.frame()

for (i in seq_along(beast_tree_files)) {
  beast_tree_file <- beast_tree_files[i]
  filename <- basename(beast_tree_file)

  # Extract clone ID from filename
  clone_match <- regmatches(filename, regexec("_([0-9]+)_tree_with_trait", filename))
  if (length(clone_match[[1]]) < 2) {
    cat("Warning: Could not extract clone ID from", filename, "\n")
    next
  }
  clone_id <- clone_match[[1]][2]

  cat("Processing clone", clone_id, paste0("(", i, "/", length(beast_tree_files), ")\n"))

  # Load BEAST tree
  beast_tree <- tryCatch({
    treeio::read.beast(beast_tree_file)
  }, error = function(e) {
    cat("  Error loading BEAST tree:", e$message, "\n")
    next
  })

  # Get corresponding true tree
  true_tree <- true_trees[[as.numeric(clone_id)]]
  if (is.null(true_tree)) {
    cat("  Warning: No true tree for clone", clone_id, "\n")
    next
  }

  # -------------------------- Tree preparation --------------------------
  # Get phylo objects
  beast_phylo <- as.phylo(beast_tree)
  true_phylo <- as.phylo(true_tree)

  # Clean tip labels
  beast_tips <- gsub("_heavy$", "", beast_phylo$tip.label)
  true_tips <- true_phylo$tip.label
  common_tips <- intersect(beast_tips, true_tips)

  # Prune trees to common tips
  beast_keep <- beast_phylo$tip.label[beast_tips %in% common_tips]
  true_keep <- true_phylo$tip.label[true_tips %in% common_tips]

  beast_pruned_phylo <- ape::keep.tip(beast_phylo, beast_keep)
  true_pruned_phylo <- ape::keep.tip(true_phylo, true_keep)

  # Clean beast labels
  beast_pruned_phylo$tip.label <- gsub("_heavy$", "", beast_pruned_phylo$tip.label)

  tip_order <- sort(common_tips)
  beast_pruned_phylo <- ape::keep.tip(beast_pruned_phylo, tip_order)
  true_pruned_phylo <- ape::keep.tip(true_pruned_phylo, tip_order)

  # -------------------------- Calculate metrics --------------------------
  # Tree heights - read BEAST height from log file
  beast_height <- as.numeric(read_tree_height_from_log(beast_tree_file))
  if (is.na(beast_height)) {
    cat("  Warning: Using phylo tree height as fallback\n")
    beast_height <- max(node.depth.edgelength(beast_pruned_phylo))
  }

  true_height <- max(node.depth.edgelength(true_pruned_phylo)) - 1

  # Tree lengths
  beast_tree_length <- sum(beast_pruned_phylo$edge.length)
  true_tree_length <- sum(true_pruned_phylo$edge.length)

  # RF distance
  rf_distance <- tryCatch({
    dowser::calcRF(di2multi(beast_pruned_phylo, tol = 5),
                   di2multi(true_pruned_phylo, tol = 5))
  }, error = function(e) {
    cat("  Warning: RF calculation failed:", e$message, "\n")
    NA
  })

  # MRCA cell type accuracy analysis
  mrca_cell_type_accuracy <- NA

  if ("location" %in% colnames(beast_tree@data) && "location" %in% colnames(true_tree@data)) {
    # Extract MRCA data using the more robust function
    beast_mrca_data <- extract_mrca_locations(beast_tree, clean_tip_label = TRUE)
    # Standardize BEAST MRCA location labels
    beast_mrca_data <- standardize_mrca_locations(beast_mrca_data)

    true_mrca_data <- extract_mrca_locations(true_tree, clean_tip_label = FALSE)

    # Filter for common tips only
    beast_mrca_filtered <- beast_mrca_data %>%
      filter(tip1 %in% common_tips & tip2 %in% common_tips)

    true_mrca_filtered <- true_mrca_data %>%
      filter(tip1 %in% common_tips & tip2 %in% common_tips)

    # Match pairs between beast and true trees
    if (nrow(beast_mrca_filtered) > 0 && nrow(true_mrca_filtered) > 0) {
      # Create matching keys
      beast_mrca_filtered$pair_key <- paste(pmin(beast_mrca_filtered$tip1, beast_mrca_filtered$tip2),
                                            pmax(beast_mrca_filtered$tip1, beast_mrca_filtered$tip2),
                                            sep = "|")
      true_mrca_filtered$pair_key <- paste(pmin(true_mrca_filtered$tip1, true_mrca_filtered$tip2),
                                           pmax(true_mrca_filtered$tip1, true_mrca_filtered$tip2),
                                           sep = "|")

      # Merge on matching pairs
      matched_data <- merge(beast_mrca_filtered, true_mrca_filtered,
                            by = "pair_key", suffixes = c("_beast", "_true"))

      if (nrow(matched_data) > 0) {
        # Calculate location accuracy
        location_matches <- ifelse(matched_data$mrca_location_beast == matched_data$mrca_location_true, 1, 0)
        mrca_cell_type_accuracy <- mean(location_matches, na.rm = TRUE)

        cat("  MRCA pairs analyzed:", nrow(matched_data), "\n")
        cat("  Location accuracy:", round(mrca_cell_type_accuracy, 3), "\n")
      }
    }
  }

  # Store results
  result <- list(
    clone_id = clone_id,
    beast_tree_height = beast_height,
    true_tree_height = true_height,
    beast_tree_length = beast_tree_length,
    true_tree_length = true_tree_length,
    rf_distance = rf_distance,
    tips_analyzed = length(common_tips),
    mrca_cell_type_accuracy = mrca_cell_type_accuracy
  )

  results_list[[clone_id]] <- result

  # Add to summary data
  summary_row <- data.frame(
    simulation_name = simulation_name,
    analysis_type = analysis_type,
    rev_suffix = rev_suffix,
    config = config_name,
    template_name = template_name,
    clone_id = clone_id,
    metric = c("beast_tree_height", "true_tree_height", "beast_tree_length", "true_tree_length",
               "rf_distance", "mrca_cell_type_accuracy", "tips_analyzed"),
    value = c(beast_height, true_height, beast_tree_length, true_tree_length,
              rf_distance, mrca_cell_type_accuracy, length(common_tips))
  )

  summary_data <- rbind(summary_data, summary_row)
  cat("  Completed clone", clone_id, "\n")
}

# -------------------------- Save results --------------------------
cat("Saving results...\n")

# Save detailed results
results_file <- file.path(config_output_dir, paste0("tree_analysis_results_", template_name, ".rds"))
saveRDS(results_list, results_file)

# Save summary CSV
summary_file <- file.path(config_output_dir, paste0("tree_analysis_summary_", template_name, ".csv"))
write.csv(summary_data, summary_file, row.names = FALSE)

# Print summary
cat("\n=== Analysis Summary ===\n")
cat("Configuration:", config_name, "\n")
cat("Template:", template_name, "\n")
cat("Clones processed:", length(results_list), "/", length(beast_tree_files), "\n")
cat("Results saved to:", results_file, "\n")
cat("Summary saved to:", summary_file, "\n")
cat("========================\n")
