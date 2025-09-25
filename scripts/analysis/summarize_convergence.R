#!/usr/bin/env Rscript
suppressMessages(library(dowser))
suppressMessages(library(airr))
suppressMessages(library(stringr))
source("scripts/utils/phylo_utilities.R")

# ------------------ Constants ------------------ #
DEFAULT_ESS_CUTOFF <- 100
DEFAULT_BURNIN <- 10
DEFAULT_NPROC <- 4
DEFAULT_NPROC_FMT <- 20
DEFAULT_MCMC_LENGTH <- 1e8

# Create clones from raw AIRR data
create_clones_from_airr <- function(config_name, simulation_name, analysis_scope, rev_suffix, model_type) {
  PROJECT_ROOT <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"

  # Find AIRR file
  raw_data_dir <- file.path(PROJECT_ROOT, simulation_name, "data/raw")
  airr_file <- file.path(raw_data_dir, config_name, "all_samples_airr.tsv")

  if (!file.exists(airr_file)) {
    cat("AIRR file not found:", airr_file, "\n")
    return(NULL)
  }

  cat("Creating clones from AIRR file:", airr_file, "\n")

  # Read and validate the input data
  data <- read_rearrangement(airr_file)
  if (!"sample_time" %in% names(data)) {
    stop("AIRR data missing 'sample_time' column")
  }

  data <- data %>%
    # Filter to heavy chains
    filter(locus == "IGH") %>%
    mutate(celltype = ifelse(celltype == "default", "GC", "other")) %>%
    # Convert sample_time to numeric
    mutate(time = as.numeric(sample_time))

  if (analysis_scope == "differentiation_analysis") {
    data$celltype <- gsub("\\s+", "", tolower(data$celltype))
  }

  # Create datasets for different analysis phases
  GC_data <- data[data$location == "germinal_center", ]
  analysis_data <- if (model_type == "gc_strict_clock") GC_data else data

  # Handle neutral vs selective evolution
  is_neutral <- grepl("_neu$", config_name)

  if (is_neutral) {
    analysis_data <- ensure_neutral_fields(analysis_data)
    GC_data       <- ensure_neutral_fields(GC_data)
  }

  # Set AIRR processing parameters based on evolution type
  airr_params <- if (is_neutral) {
    list(filterstop = FALSE, use_regions = FALSE, germ = "germline_alignment")
  } else {
    list(filterstop = TRUE,  use_regions = TRUE,  germ = "germline_alignment")
  }

  cat("Processing data: Evolution =", if (is_neutral) "neutral" else "selective", "\n")

  # Format clones with appropriate parameters
  clones <- formatClones(
    analysis_data,
    traits      = c("location", "time", "celltype"),
    germ        = "germline_alignment",
    chain       = "H",
    nproc       = DEFAULT_NPROC_FMT,
    collapse    = TRUE,
    minseq      = 3,
    filterstop  = airr_params$filterstop,
    use_regions = airr_params$use_regions,
    germ        = airr_params$germ,
    v_call      = "v_call",
    j_call      = "j_call",
    junc_len    = "junction_length"
  )

  # Save processed clones for future use
  processed_data_dir <- file.path(PROJECT_ROOT, simulation_name, "data/processed", analysis_scope)
  dir.create(processed_data_dir, recursive = TRUE, showWarnings = FALSE)
  clones_file <- file.path(processed_data_dir, paste0(config_name, "_", rev_suffix, "_filtered_clones.rds"))
  saveRDS(clones, clones_file)
  cat("Saved recreated clones to:", clones_file, "\n")

  return(clones)
}

# ------------------ Main Script ------------------ #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_convergence_summary.R <simulation_name> <rev_suffix> <analysis_type>")
}

simulation_name <- args[1]
rev_suffix <- args[2]
analysis_type <- args[3]

# Setup paths
PROJECT_ROOT <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
BEAST_DIR <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast/bin"

# Read job list CSV
configs_dir <- file.path(PROJECT_ROOT, simulation_name, "configs")
job_list_file <- file.path(configs_dir, paste0("tree_analysis_jobs_", analysis_type, "_", rev_suffix, ".csv"))

if (!file.exists(job_list_file)) {
  stop("Job list CSV not found: ", job_list_file)
}

cat("=== Reading job list from:", job_list_file, "===\n")
jobs <- read.csv(job_list_file, stringsAsFactors = FALSE)
cat("Found", nrow(jobs), "jobs\n")

# Setup directories
processed_data_dir <- file.path(PROJECT_ROOT, simulation_name, "data/processed", analysis_type)
dowser_processed_trees_dir <- file.path(processed_data_dir, "dowser_processed_trees", rev_suffix)
results_base_dir <- file.path(PROJECT_ROOT, simulation_name, "results", analysis_type, rev_suffix)

# Create summary data frame
all_summaries <- data.frame()

# Process each unique config-template combination
unique_jobs <- unique(jobs[, c("config_name", "template_name", "model_type")])

cat("\n=== Processing", nrow(unique_jobs), "unique config-template combinations ===\n")

for (i in seq_len(nrow(unique_jobs))) {
  config_name <- unique_jobs$config_name[i]
  template_name <- unique_jobs$template_name[i]
  model_type <- unique_jobs$model_type[i]

  cat("\n--- Processing:", config_name, "-", template_name, "(", model_type, ") ---\n")

  # Try to load existing dowser files
  trees <- NULL
  trees_source <- NULL

  # Try standard naming: config_name_template_name_all_trees.rds
  if (is.null(trees)) {
    standard_file <- file.path(dowser_processed_trees_dir, paste0(config_name, "_", template_name, "_all_trees.rds"))
    if (file.exists(standard_file)) {
      cat("Found standard dowser file:", standard_file, "\n")
      trees <- readRDS(standard_file)
      trees_source <- "standard_dowser"
    }
  }

  # Check if loaded trees have convergence info
  if (!is.null(trees) && is.data.frame(trees) && "below_ESS" %in% names(trees)) {
    cat("Loaded trees with existing convergence information from", trees_source, "\n")
    cat("Trees data frame has", nrow(trees), "rows and columns:", paste(names(trees), collapse = ", "), "\n")
  } else {
    # Fallback to readBEAST if no dowser files or missing convergence info
    if (is.null(trees)) {
      cat("No dowser files found, using readBEAST fallback...\n")
    } else {
      cat("Dowser file missing convergence info, using readBEAST fallback...\n")
    }

    # Check for clones file first
    clones_file <- file.path(processed_data_dir, paste0(config_name, "_", rev_suffix, "_filtered_clones.rds"))
    clones_data <- NULL

    if (file.exists(clones_file)) {
      cat("Found existing clones file:", clones_file, "\n")
      clones_data <- readRDS(clones_file)
    } else {
      cat("Clones file missing, recreating from raw AIRR data...\n")
      clones_data <- create_clones_from_airr(config_name, simulation_name, analysis_type, rev_suffix, model_type)

      if (is.null(clones_data)) {
        cat("Failed to create clones for", config_name, "\n")
        next
      }
    }

    # Find BEAST output directory
    beast_output_dir <- file.path(PROJECT_ROOT, simulation_name, "results", analysis_type, rev_suffix, model_type, "beast_raw_output", config_name)

    if (!dir.exists(beast_output_dir)) {
      cat("Warning: BEAST output directory not found:", beast_output_dir, "\n")
      next
    }

    # Call readBEAST
    tryCatch({
      trees <- readBEAST(
        clones = clones_data,
        dir = beast_output_dir,
        id = template_name,
        beast = BEAST_DIR,
        burnin = DEFAULT_BURNIN,
        nproc = DEFAULT_NPROC,
        quiet = 0,
        full_posterior = FALSE,
        asr = FALSE,
        low_ram = TRUE
      )

      cat("Created", nrow(trees), "tree objects using readBEAST\n")
      trees_source <- "readBEAST"

      # Compute below_ESS for readBEAST output
      cat("Computing ESS convergence for readBEAST trees...\n")
      ignore_params <- build_ignore_vector(template_name, model_type, analysis_type)

      trees$below_ESS <- integer(nrow(trees))
      trees$failing_params <- character(nrow(trees))

      for (j in seq_len(nrow(trees))) {
        params <- trees@info$parameters[[j]]

        if (is.null(params) || !"ESS" %in% names(params)) {
          trees$below_ESS[j] <- 1
          trees$failing_params[j] <- "no_ESS_data"
          next
        }

        # Filter out ignored parameters
        for (ignore_param in ignore_params) {
          params <- params[!grepl(ignore_param, params$item), ]
        }

        # Check ESS values
        below_ess_idx <- which(params$ESS < DEFAULT_ESS_CUTOFF)
        trees$below_ESS[j] <- length(below_ess_idx)
        trees$failing_params[j] <- if (length(below_ess_idx) > 0) {
          paste(params$item[below_ess_idx], collapse = ";")
        } else {
          ""
        }
      }

      # Save the readBEAST results with convergence info
      dir.create(dowser_processed_trees_dir, recursive = TRUE, showWarnings = FALSE)
      readbeast_file <- file.path(dowser_processed_trees_dir, paste0(config_name, "_", template_name, "_all_trees.rds"))
      saveRDS(trees, readbeast_file)
      cat("Saved readBEAST trees with ESS info to:", readbeast_file, "\n")

    }, error = function(e) {
      cat("Error calling readBEAST:", e$message, "\n")
      next
    })
  }

  if (is.null(trees) || !"below_ESS" %in% names(trees)) {
    cat("Warning: Could not load or create trees with convergence info for", config_name, template_name, "\n")
    next
  }

  # ------------------ Convergence analysis and filtering ------------------ #
  # Validate trees output
  if (!"below_ESS" %in% names(trees)) {
    stop("Expected column `below_ESS` not found in trees output")
  }

  # Identify unconverged clones
  unconverged_idx <- which(trees$below_ESS > 0)
  unconverged_clone_ids <- trees$clone_id[unconverged_idx]

  # Create converged trees dataset
  if (length(unconverged_idx) > 0) {
    trees_converged <- trees[setdiff(seq_len(nrow(trees)), unconverged_idx), ]
    cat("Filtered to", nrow(trees_converged), "converged clones out of", nrow(trees), "total\n")
  } else {
    trees_converged <- trees
    cat("All", nrow(trees), "clones converged\n")
  }

  # Log convergence summary
  cat(sprintf("Final: %d/%d clones converged (ESS ≥ %d for all non-ignored parameters)\n",
              nrow(trees_converged), nrow(trees), DEFAULT_ESS_CUTOFF))

  # ------------------ Results saving ------------------ #
  cat("=== Saving Analysis Results ===\n")

  # Save converged trees only
  converged_trees_file <- file.path(dowser_processed_trees_dir, paste0(config_name, "_", template_name, "_converged_trees.rds"))
  saveRDS(trees_converged, converged_trees_file)
  cat("Converged trees saved to:", converged_trees_file, "\n")

  # Create detailed convergence summary
  ignore_params <- build_ignore_vector(template_name, model_type, analysis_type)
  convergence_summary <- data.frame(
    config_name = config_name,
    template_id = template_name,
    model_type = model_type,
    analysis_scope = analysis_type,
    reversible = (rev_suffix == "rev"),
    total_clones = nrow(trees),
    converged_clones = nrow(trees_converged),
    unconverged_clones = length(unconverged_idx),
    unconverged_clone_ids = if (length(unconverged_idx) > 0) paste(unconverged_clone_ids, collapse = ";") else "",
    final_iteration = 10,
    mcmc_length = DEFAULT_MCMC_LENGTH,
    ess_cutoff = DEFAULT_ESS_CUTOFF,
    ignore_params = paste(ignore_params, collapse = ";"),
    trees_source = trees_source,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )

  all_summaries <- rbind(all_summaries, convergence_summary)
}

# Save convergence summaries by model type
cat("\n=== Saving convergence summaries ===\n")

for (model_type in unique(all_summaries$model_type)) {
  model_summaries <- all_summaries[all_summaries$model_type == model_type, ]
  summary_dir <- file.path(results_base_dir, model_type, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  # Update combined convergence summary
  combined_summary_file <- file.path(summary_dir, paste0(model_type, "_convergence_summary.csv"))

  if (file.exists(combined_summary_file)) {
    # Read existing summary and update
    existing_summary <- read.csv(combined_summary_file, stringsAsFactors = FALSE)

    # Remove any existing entries for configs/templates we just processed
    for (i in seq_len(nrow(model_summaries))) {
      existing_summary <- existing_summary[!(existing_summary$config_name == model_summaries$config_name[i] &
                                               existing_summary$template_id == model_summaries$template_id[i]), ]
    }

    # Add new entries
    updated_summary <- rbind(existing_summary, model_summaries)
    write.csv(updated_summary, combined_summary_file, row.names = FALSE)
    cat("Updated combined convergence summary:", combined_summary_file, "\n")
  } else {
    # Create new combined summary
    write.csv(model_summaries, combined_summary_file, row.names = FALSE)
    cat("Created new combined convergence summary:", combined_summary_file, "\n")
  }
}

cat("\n=== Final Summary ===\n")
cat("Total configurations processed:", length(unique(all_summaries$config_name)), "\n")
cat("Total templates processed:", length(unique(all_summaries$template_id)), "\n")
cat("Total converged clones:", sum(all_summaries$converged_clones), "\n")
cat("Total unconverged clones:", sum(all_summaries$unconverged_clones), "\n")

# Show source breakdown
if ("trees_source" %in% names(all_summaries)) {
  source_summary <- table(all_summaries$trees_source)
  cat("Trees source breakdown:\n")
  for (src in names(source_summary)) {
    cat("  ", src, ":", source_summary[src], "configs\n")
  }
}

# Print configs with failed clones
failed_configs <- unique(all_summaries$config_name[all_summaries$unconverged_clones > 0])
if (length(failed_configs) > 0) {
  cat("Configs with failed clones:", paste(failed_configs, collapse = ", "), "\n")
} else {
  cat("All configs converged successfully\n")
}

cat("=== Convergence summary generation completed ===\n")