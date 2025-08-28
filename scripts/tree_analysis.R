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
if (length(args) < 3) {
  cat("Usage: Rscript tree_analysis.R job_list_csv array_task_id output_dir\n")
  quit(status = 1)
}

job_list_file <- args[1]
array_task_id <- as.numeric(args[2])
output_dir <- args[3]

# Read CSV and get the specific job
job_data <- read.csv(job_list_file)
job_row <- job_data[array_task_id, ]

config_name <- job_row$config_name
template_name <- job_row$template_name
true_tree_file <- job_row$true_tree_file
beast_files <- strsplit(job_row$beast_files, ";")[[1]]

cat("=== Tree Analysis Started ===\n")
cat("Config:", config_name, "\n")
cat("Template:", template_name, "\n")
cat("BEAST files:", length(beast_files), "\n")

# Create output directory
config_output_dir <- file.path(output_dir, "results", config_name)
dir.create(config_output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to read TreeHeight from BEAST log file
read_tree_height_from_log <- function(beast_file) {
  # Construct log file path
  log_tsv_file <- gsub("_tree_with_trait\\.tree$", "_log.tsv", beast_file)

  if (!file.exists(log_tsv_file)) {
    cat("  Warning: Log file not found:", basename(log_tsv_file), "\n")
    return(NA)
  }
  # log <- fread(file_name, skip = 3, header=TRUE) %>% as.data.frame()
  # Read TSV file
  lines <- readLines(log_tsv_file, warn = FALSE)
  tree_height_line <- grep("^TreeHeight\\s+", lines, value = TRUE)
  tree_height_mean <- unlist(strsplit(tree_height_line[1], "\\s+"))[2]

  return(tree_height_mean)
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

for (i in seq_along(beast_files)) {
  beast_file <- beast_files[i]
  filename <- basename(beast_file)

  # Extract clone ID from filename
  clone_match <- regmatches(filename, regexec("_([0-9]+)_tree_with_trait", filename))
  if (length(clone_match[[1]]) < 2) {
    cat("Warning: Could not extract clone ID from", filename, "\n")
    next
  }
  clone_id <- clone_match[[1]][2]

  cat("Processing clone", clone_id, paste0("(", i, "/", length(beast_files), ")\n"))

  # Load BEAST tree
  beast_tree <- tryCatch({
    treeio::read.beast(beast_file)
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
  beast_height <- as.numeric(read_tree_height_from_log(beast_file))
  if (is.na(beast_height)) {
    cat("  Warning: Using phylo tree height as fallback\n")
    beast_height <- max(node.depth.edgelength(beast_pruned_phylo))
  }

  true_height <- max(node.depth.edgelength(true_pruned_phylo)) - 1

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
  mrca_height_correlation <- NA

  if ("location" %in% colnames(beast_tree@data) && "location" %in% colnames(true_tree@data)) {
    # Get node heights
    beast_node_heights <- node.depth.edgelength(beast_pruned_phylo)
    true_node_heights <- node.depth.edgelength(true_pruned_phylo)

    # MRCA comparison for paired tips
    n_tips <- length(common_tips)
    if (n_tips >= 2) {
      mrca_beast_heights <- numeric()
      mrca_true_heights <- numeric()
      location_matches <- numeric()

      for (i in 1:(n_tips - 1)) {
        for (j in (i + 1):n_tips) {
          # Get MRCA nodes
          beast_mrca <- getMRCA(beast_pruned_phylo, c(common_tips[i], common_tips[j]))
          true_mrca <- getMRCA(true_pruned_phylo, c(common_tips[i], common_tips[j]))

          if (!is.null(beast_mrca) && !is.null(true_mrca)) {
            mrca_beast_heights <- c(mrca_beast_heights, beast_node_heights[beast_mrca])
            mrca_true_heights <- c(mrca_true_heights, true_node_heights[true_mrca])

            # Check location match
            beast_loc <- beast_tree@data[as.numeric(beast_tree@data$node) == beast_mrca, "location"]
            true_loc <- true_tree@data[as.numeric(true_tree@data$node) == true_mrca, "location"]

            if (length(beast_loc) > 0 && length(true_loc) > 0) {
              location_matches <- c(location_matches,
                                    ifelse(beast_loc[1] == true_loc[1], 1, 0))
            }
          }
        }
      }

      # Calculate metrics
      if (length(mrca_beast_heights) > 0) {
        mrca_cell_type_accuracy <- mean(location_matches, na.rm = TRUE)
        mrca_height_correlation <- cor(mrca_beast_heights, mrca_true_heights, use = "complete.obs")
      }
    }
  }

  # Store results
  result <- list(
    clone_id = clone_id,
    beast_tree_height = beast_height,
    true_tree_height = true_height,
    rf_distance = rf_distance,
    tips_analyzed = length(common_tips),
    mrca_cell_type_accuracy = mrca_cell_type_accuracy,
    mrca_height_correlation = mrca_height_correlation
  )

  results_list[[clone_id]] <- result

  # Add to summary data
  summary_row <- data.frame(
    config = config_name,
    template_name = template_name,
    clone_id = clone_id,
    metric = c("beast_tree_height", "true_tree_height", "rf_distance",
               "mrca_cell_type_accuracy", "mrca_height_correlation", "tips_analyzed"),
    value = c(beast_height, true_height, rf_distance,
              mrca_cell_type_accuracy, mrca_height_correlation, length(common_tips))
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
cat("Clones processed:", length(results_list), "/", length(beast_files), "\n")
cat("Results saved to:", results_file, "\n")
cat("Summary saved to:", summary_file, "\n")
cat("========================\n")
