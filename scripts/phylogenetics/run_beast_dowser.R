# devtools::install_github("immcantation/dowser@bug_fixes")
suppressMessages(library(dowser))
suppressMessages(library(airr))
suppressMessages(library(stringr))
source("scripts/analysis/tree_functions.R")

# ------------------ Constants and configuration ------------------ #
# Default parameters for BEAST and dowser analysis
DEFAULT_ESS_CUTOFF   <- 200
DEFAULT_MAX_ITER     <- 10
DEFAULT_PERMUTATIONS <- 100000
DEFAULT_MCMC_LENGTH  <- 1e8
DEFAULT_LOG_TARGET   <- 2000
DEFAULT_NPROC_FMT    <- 20
DEFAULT_NPROC_BEAST  <- 20
DEFAULT_SEED         <- 12345

PROJECT_ROOT <- "/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
SCRATCH_ROOT <- "/dartfs-hpc/scratch/f0070d5"

# Template mapping for different models
TEMPLATE_MAP <- c(
  # GC strict clock model
  "StrictClock_Standard"                                   = "StrictClock/StrictClock_Standard_EmpFreq.xml",
  # TyCHE models
  "ExpectedOccupancy_FixedClockRates"                      = "TypeLinked/TraitLinkedExpectedOccupancy_FixedTraitClockRates_EmpFreq.xml",
  "ExpectedOccupancy_EstClockRates"                        = "TypeLinked/TraitLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
  "MixedSwitch_FixedClockRates"                            = "TypeLinked/TraitLinkedMixedSwitch_FixedTraitClockRates_EmpFreq.xml",
  "MixedSwitch_EstClockRates"                              = "TypeLinked/TraitLinkedMixedSwitch_EstTraitClockRates_EmpFreq.xml",
  "InstantSwitch_EstClockRates"                            = "TypeLinked/TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq.xml",
  # Competing models
  "StrictClock_AncestralReconstruction"                    = "StrictClock/StrictClock_AncestralReconstruction_EmpFreq.xml",
  "UCRelaxedClock_AncestralReconstruction"                 = "UCLD/UCRelaxedClock_AncestralReconstruction_EmpFreq.xml",
  # TyCHE 3-state model for differentiation analysis
  "TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq_3state" = "TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq_3state.xml"
)

# ------------------ Args ------------------ #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 9) {
  task_id <- as.numeric(args[1])
  combinations_file <- args[2]
  simulation_run <- args[3]
  analysis_scope <- args[4]
  model_type <- args[5]
  reversible <- as.logical(args[6])
  beast_dir <- args[7]
  nproc_arg <- as.numeric(args[8])
  seed <- as.numeric(args[9])
} else {
  stop("Usage: Rscript run_beast_dowser.R <task_id> <combinations_file> <simulation_run> <analysis_scope> <model_type> <reversible> <beast_dir> <nproc> <seed>")
}

cat("=== BEAST Dowser Job", task_id, "Started at", as.character(Sys.time()), "===\n")

# Read job combinations
combinations <- read.csv(combinations_file, stringsAsFactors = FALSE)
job_params <- combinations[task_id, ]

# Extract job parameters
config_name <- job_params$config_name
template_id <- job_params$template_id
include_germline <- job_params$germline
airr_file <- job_params$airr_file

# Handle time subset for sub_analysis
time_filter <- NULL
if (analysis_scope == "sub_analysis" && "time_subset" %in% names(job_params)) {
  time_subset_str <- job_params$time_subset
  time_filter <- as.numeric(unlist(strsplit(time_subset_str, "_")))
  cat("Time subset for sub_analysis:", paste(time_filter, collapse = ", "), "\n")
}

# Configure analysis-specific parameters
if (analysis_scope == "main_analysis") {
  trait_param <- "location"
  template_dir_suffix <- "custom"
} else if (analysis_scope == "sub_analysis") {
  trait_param <- "location"
  template_dir_suffix <- "custom"
} else if (analysis_scope == "differentiation_analysis") {
  trait_param <- "celltype"
  template_dir_suffix <- if (model_type == "gc_strict_clock") "custom" else "ken_templates"
}

# ------------------ Helpers ------------------ #
# Create directory if it doesn't exist
create_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Log key-value pairs in a formatted way
logKVs <- function(title, kvs) {
  cat(strrep("=", 15), title, strrep("=", 15), "\n", sep = " ")
  for (nm in names(kvs))
  cat(sprintf("%-12s: %s\n", nm, kvs[[nm]]))
}

logKVs("Job", c(
  Config           = config_name,
  Template         = template_id,
  Reversible       = reversible,
  Step             = model_type,
  NPROC            = nproc_arg
))

# Template path setup
template_dir <- file.path(PROJECT_ROOT, "xml-writer/templates", template_dir_suffix)
template_file <- TEMPLATE_MAP[[template_id]]
if (is.null(template_file)) {
  stop("Unknown template ID: ", template_id)
}
template_path <- file.path(template_dir, template_file)
if (!file.exists(template_path)) {
  stop("Template file not found: ", template_path)
}

# ------------------ Paths ------------------ #
# Determine reversibility suffix
rev_suffix <- if (reversible) "rev" else "irrev"

# Setup directory paths consistent with bash script
base_data_dir <- file.path(PROJECT_ROOT, simulation_run, "data")
raw_data_dir <- file.path(base_data_dir, "raw")
processed_data_dir <- file.path(base_data_dir, "processed")
analysis_processed_data_dir <- create_dir(file.path(processed_data_dir, analysis_scope))
dowser_processed_trees_dir <- create_dir(file.path(analysis_processed_data_dir, "dowser_processed_trees", rev_suffix))

results_base_dir <- file.path(PROJECT_ROOT, simulation_run, "results", analysis_scope, rev_suffix, model_type)
gc_strict_results_base_dir <- file.path(PROJECT_ROOT, simulation_run, "results", analysis_scope, rev_suffix, "gc_strict_clock")

# Output directories
beast_raw_output_dir <- create_dir(file.path(SCRATCH_ROOT, simulation_run, analysis_scope, rev_suffix, model_type, config_name))
summary_dir <- create_dir(file.path(results_base_dir, "summary"))

# Model-specific directories
if (model_type == "gc_strict_clock") {
  correlation_dir <- create_dir(file.path(results_base_dir, "correlation_analysis"))
  clock_rates_dir <- create_dir(file.path(results_base_dir, "clock_rates"))
}

# Analysis-specific directories
if (analysis_scope == "sub_analysis" && !is.null(time_filter)) {
  time_suffix <- paste(time_filter, collapse = "_")

  # Create time-specific subdirectory under dowser_processed_trees and results directory
  time_specific_trees_dir <- create_dir(file.path(dowser_processed_trees_dir, paste0("time_", time_suffix)))
  beast_raw_output_dir <- create_dir(file.path(SCRATCH_ROOT, simulation_run, analysis_scope,
                                               rev_suffix, model_type, paste0("time_", time_suffix), config_name))
}

logKVs("Directory Paths", c(
  Raw_data         = raw_data_dir,
  Processed_data   = processed_data_dir,
  Analysis_Processed_Data = analysis_processed_data_dir,
  Dowser_trees     = dowser_processed_trees_dir,
  Beast_output     = beast_raw_output_dir,
  Results_base     = results_base_dir,
  Template_path    = template_path
))

# ------------------ Data loading and preprocessing ------------------ #
# Read and validate the input data
if (!file.exists(airr_file)) {
  stop("AIRR file not found: ", airr_file)
}

data <- read_rearrangement(airr_file)
if (!"sample_time" %in% names(data)) {
  stop("AIRR data missing 'sample_time' column")
}

# Convert sample_time to numeric
data$time <- as.numeric(as.character(data$sample_time))

# Filter to heavy chains and specific clones for analysis
data <- data[data$locus == "IGH", ]
if (model_type == "gc_strict_clock") {
  cat("Running clone IDs 1-20 for germinal center strict clock analysis\n")
} else {
  cat("Running clone IDs 1-20 for", model_type, "analysis\n")
}

if (analysis_scope == "differentiation_analysis") {
  data$celltype <- gsub("\\s+", "", tolower(data$celltype))
}
# Apply time filter if specified for sub_analysis
if (!is.null(time_filter)) {
  seq_count <- nrow(data)
  data <- data[data$time %in% time_filter, ]
  cat("Time filtering: ", seq_count, "→", nrow(data), "sequences (times:", paste(time_filter, collapse = ", "), ")\n")
}

# Create datasets for different analysis phases
GC_data <- data[data$location == "germinal_center", ]
analysis_data <- if (model_type == "gc_strict_clock") GC_data else data

# Handle neutral vs selective evolution
is_neutral <- grepl("_neu$", config_name)
is_selective <- grepl("_sel$", config_name)

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

logKVs("AIRR Data Processing Summary", c(
  Evolution   = if (is_neutral) "neutral" else "selective",
  Filter_Stop     = airr_params$filterstop,
  Use_Regions     = airr_params$use_regions,
  Germline_Field  = airr_params$germ
))

# ------------------ Clone formatting and tree building ------------------ #
# Format clones with appropriate parameters
clones <- formatClones(
  analysis_data,
  traits      = c("location", "time", "celltype"),
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

# Save processed clones
clones_file <- file.path(analysis_processed_data_dir, paste0(config_name, "_", rev_suffix, "_filtered_clones.rds"))
if (analysis_scope == "sub_analysis" && !is.null(time_filter)) {
  clones_file <- file.path(analysis_processed_data_dir, paste0("time_", time_suffix, "_", config_name, "_filtered_clones.rds"))
}
saveRDS(clones, clones_file)
cat("Saved clones to:", clones_file, "\n")

# Perform correlation analysis for GC strict clock models
if (model_type == "gc_strict_clock") {
  cat("=== Running Correlation Analysis ===\n")

  trees_for_correlation <- getTrees(clones, build = "pml", sub_model = "HKY")
  cat("Built", nrow(trees_for_correlation), "trees for correlation analysis\n")

  correlation_results <- correlationTest(trees_for_correlation, permutations = DEFAULT_PERMUTATIONS)
  correlation_tbl <- correlation_results[, c("clone_id", "locus", "slope", "p")]

  slopes <- correlation_tbl$slope

  correlation_file <- file.path(correlation_dir, paste0(config_name, "_correlation_results.csv"))
  write.csv(correlation_tbl, correlation_file, row.names = FALSE)
  cat("Correlation results saved to:", correlation_file, "\n")

  cat("=== Correlation Analysis Summary ===\n")
  cat("Total clones analyzed:", nrow(correlation_tbl), "\n")
  if (length(slopes)) {
    cat("Correlation slopes mean/median:", round(mean(slopes), 6), "/", round(median(slopes), 6), "\n")
  } else {
    cat("WARNING: No correlation slopes found\n")
  }
  cat("===================================\n")
}

# ------------------ Model-specific parameter setup ------------------ #
RATE_INDICATORS <- if (reversible) "1 1" else "1 0"
KAPPA_PRIOR_M   <- "0.67"
KAPPA_PRIOR_S   <- "0.2"

# TyCHE model parameters
TRAIT_RATE_MEAN_1 <- NULL
TRAIT_RATE_MEAN_2 <- NULL
TRAIT_RATE_MEAN_3 <- NULL
TRAIT_RATE_SIGMA_1 <- NULL
TRAIT_RATE_SIGMA_2 <- NULL
TRAIT_RATE_SIGMA_3 <- NULL

# Strict clock and UCLD model parameters
CLOCK_RATE_INIT   <- NULL

# UCLD model parameters
UCLD_SIGMA_INIT   <- NULL

# Transition rate parameters
TRANSITION_RATE_ALPHA_1 <- NULL
TRANSITION_RATE_BETA_1 <- NULL
TRANSITION_RATE_ALPHA_2 <- NULL
TRANSITION_RATE_BETA_2 <- NULL
TRANSITION_RATE_ALPHA <- NULL
TRANSITION_RATE_BETA <- NULL

# Differentiation analysis-specific parameters
GEO_RATE_ALPHA <- NULL
GEO_RATE_BETA <- NULL
INITIAL_STATE <- NULL

# Extract config-specific mean GC clock rates for TyCHE and competing models during all analysis scopes
# and GC stict clock model for differentiation timing analysis
if (model_type %in% c("tyche_models", "competing_models")) {
  # Determine which analysis scope to get GC rates from
  if (analysis_scope == "differentiation_analysis") {
    # Use main_analysis GC rates for differentiation analysis
    gc_rates_analysis_scope <- "main_analysis"
    cat("Using GC clock rates from main_analysis for differentiation_analysis\n")
  } else {
    # Use own analysis scope GC rates (for sub_analysis)
    gc_rates_analysis_scope <- analysis_scope
    cat("Using GC clock rates from", analysis_scope, "\n")
  }

  # Build path to appropriate GC results
  gc_results_dir <- file.path(PROJECT_ROOT, simulation_run, "results", gc_rates_analysis_scope, rev_suffix, "gc_strict_clock")
  clock_rates_dir <- file.path(gc_results_dir, "clock_rates")
  individual_mean_rates_file <- file.path(clock_rates_dir, paste0(config_name, "_individual_mean_gc_clock_rates.csv"))

  if (file.exists(individual_mean_rates_file)) {
    mean_rates_data <- read.csv(individual_mean_rates_file)
    matching_config <- mean_rates_data[mean_rates_data$config_name == config_name, ]
    if (nrow(matching_config) > 0) {
      mean_gc_clock_rate <- matching_config$mean_clock_rate
    } else {
      stop("No matching GC strict clock config found for: ", config_name)
    }
  } else {
    stop("GC strict clock mean rates CSV not found: ", individual_mean_rates_file)
  }
}

# Configure parameters based on analysis scope and model type
if (analysis_scope %in% c("main_analysis", "sub_analysis")) {
  if (model_type == "gc_strict_clock") {
    # Simple strict clock on GC data
    trait_param <- NULL
    CLOCK_RATE_INIT  <- "0.001"

  } else if (model_type == "tyche_models") {
    # TyCHE models on all B cell populations
    TRAIT_RATE_MEAN_1 <- as.character(mean_gc_clock_rate)
    TRAIT_RATE_MEAN_2 <- "0.000001"

    # Transition rate parameters for TyCHE
    TRANSITION_RATE_ALPHA_1 <- "0.1"
    TRANSITION_RATE_BETA_1 <- "1.0"
    TRANSITION_RATE_ALPHA_2 <- "0.1"
    TRANSITION_RATE_BETA_2 <- "1.0"

  } else if (model_type == "competing_models") {
    # Competing models (Strict + UCLD) on all B cell populations
    CLOCK_RATE_INIT <- as.character(mean_gc_clock_rate)

    # Transition rate parameters for competing models
    TRANSITION_RATE_ALPHA <- "0.1"
    TRANSITION_RATE_BETA <- "1.0"

    # UCLD-specific parameters
    if (grepl("UCRelaxedClock", template_id)) {
      UCLD_SIGMA_INIT <- "0.5"
    }
  }

  # Calculate sigma parameters for trait rates (if applicable)
  # Default priors inlined (strong_GC_weak_mbc)
  if (!is.null(TRAIT_RATE_MEAN_1)) {
    TRAIT_RATE_SIGMA_1 <- as.character(0.01 * as.numeric(TRAIT_RATE_MEAN_1))
    TRAIT_RATE_SIGMA_2 <- "0.001"
  }

} else if (analysis_scope == "differentiation_analysis") {
  if (model_type == "gc_strict_clock") {
    # Simple strict clock on GC data
    trait_param <- NULL
    CLOCK_RATE_INIT  <- "0.001"

  } else if (grepl("3state", template_id)) {
    # 3-state TyCHE model parameters
    TRAIT_RATE_MEAN_1 <- as.character(mean_gc_clock_rate)
    TRAIT_RATE_MEAN_2 <- format(1e-6, scientific = F)
    TRAIT_RATE_MEAN_3 <- format(1e-6, scientific = F)
    TRAIT_RATE_SIGMA_1 <- as.character(0.01 * as.numeric(TRAIT_RATE_MEAN_1))
    TRAIT_RATE_SIGMA_2 <- "0.001"
    TRAIT_RATE_SIGMA_3 <- "0.001"
    RATE_INDICATORS <- "1 1 0 0 0 0"
    TRANSITION_RATE_ALPHA <- "0.1"
    TRANSITION_RATE_BETA <- "1.0"
    GEO_RATE_ALPHA <- "0.01"
    GEO_RATE_BETA <- "1.0"
    INITIAL_STATE <- "0"
  }
}

logKVs("BEAST Parameters Summary", c(
  Model_Type          = model_type,
  Trait_Parameter     = if (is.null(trait_param)) "NULL" else trait_param,
  TRAIT_RATE_MEAN_1   = if (is.null(TRAIT_RATE_MEAN_1)) "NULL" else TRAIT_RATE_MEAN_1,
  TRAIT_RATE_MEAN_2   = if (is.null(TRAIT_RATE_MEAN_2)) "NULL" else TRAIT_RATE_MEAN_2,
  TRAIT_RATE_MEAN_3   = if (is.null(TRAIT_RATE_MEAN_3)) "NULL" else TRAIT_RATE_MEAN_3,
  CLOCK_RATE_INIT     = if (is.null(CLOCK_RATE_INIT)) "NULL" else CLOCK_RATE_INIT,
  UCLD_SIGMA_INIT     = if (is.null(UCLD_SIGMA_INIT)) "NULL" else UCLD_SIGMA_INIT,
  RATE_INDICATORS     = RATE_INDICATORS,
  Template_ID         = template_id
))

# ------------------ Run getTimeTreesIterate ------------------ #
# Create XML files and build time trees with iterative improvement

# Setup convergence parameters
ess_cutoff <- DEFAULT_ESS_CUTOFF
max_iter   <- DEFAULT_MAX_ITER
ignore_params <- build_ignore_vector(template_id, model_type, analysis_scope)

cat("=== Running getTimeTreesIterate ===\n")
cat("ESS cutoff:", ess_cutoff, "\n")
cat("Max iterations:", max_iter, "\n")
cat("Ignore parameters:", paste(ignore_params, collapse = ", "), "\n")

# Run iterative BEAST analysis
trees <- getTimeTreesIterate(
  clones             = clones,
  iterations         = max_iter,
  ess_cutoff         = ess_cutoff,
  ignore             = ignore_params,
  quiet              = 0,
  # getTimeTrees()
  beast              = beast_dir,
  trait              = trait_param,
  time               = "time",
  dir                = beast_raw_output_dir,
  id                 = template_id,
  template           = template_path,
  nproc              = nproc_arg,
  include_germline   = include_germline,
  mcmc_length        = DEFAULT_MCMC_LENGTH,
  log_every          = "auto",
  log_target         = DEFAULT_LOG_TARGET,
  seed               = seed + (task_id * 10000),
  # Model-specific parameters
  TRAIT_RATE_MEAN_1  = TRAIT_RATE_MEAN_1,
  TRAIT_RATE_MEAN_2  = TRAIT_RATE_MEAN_2,
  TRAIT_RATE_SIGMA_1 = TRAIT_RATE_SIGMA_1,
  TRAIT_RATE_SIGMA_2 = TRAIT_RATE_SIGMA_2,
  RATE_INDICATORS    = RATE_INDICATORS,
  CLOCK_RATE_INIT    = CLOCK_RATE_INIT,
  KAPPA_PRIOR_M      = KAPPA_PRIOR_M,
  KAPPA_PRIOR_S      = KAPPA_PRIOR_S,
  UCLD_SIGMA_INIT    = UCLD_SIGMA_INIT,
  # Transition rate parameters
  TRANSITION_RATE_ALPHA_1 = TRANSITION_RATE_ALPHA_1,
  TRANSITION_RATE_BETA_1  = TRANSITION_RATE_BETA_1,
  TRANSITION_RATE_ALPHA_2 = TRANSITION_RATE_ALPHA_2,
  TRANSITION_RATE_BETA_2  = TRANSITION_RATE_BETA_2,
  TRANSITION_RATE_ALPHA   = TRANSITION_RATE_ALPHA,
  TRANSITION_RATE_BETA    = TRANSITION_RATE_BETA,
  # 3-state TyCHE-specific parameters
  TRAIT_RATE_MEAN_3  = TRAIT_RATE_MEAN_3,
  TRAIT_RATE_SIGMA_3 = TRAIT_RATE_SIGMA_3,
  GEO_RATE_ALPHA     = GEO_RATE_ALPHA,
  GEO_RATE_BETA      = GEO_RATE_BETA,
  INITIAL_STATE      = INITIAL_STATE
)

cat("=== getTimeTreesIterate completed ===\n")

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
            nrow(trees_converged), nrow(trees), ess_cutoff))

# ------------------ Results saving ------------------ #
cat("=== Saving Analysis Results ===\n")
if (analysis_scope == "sub_analysis" && !is.null(time_filter)) {
  dowser_processed_trees_dir <- time_specific_trees_dir
  dowser_processed_trees_dir <- time_specific_trees_dir
}

# Save all trees (including unconverged)
all_trees_file <- file.path(dowser_processed_trees_dir, paste0(config_name, "_", template_id, "_all_trees.rds"))
saveRDS(trees, all_trees_file)
cat("All trees saved to:", all_trees_file, "\n")

# Save converged trees only
converged_trees_file <- file.path(dowser_processed_trees_dir, paste0(config_name, "_", template_id, "_converged_trees.rds"))
saveRDS(trees_converged, converged_trees_file)
cat("Converged trees saved to:", converged_trees_file, "\n")

# Create detailed convergence summary
convergence_summary <- data.frame(
  config_name        = config_name,
  template_id        = template_id,
  model_type         = model_type,
  analysis_scope     = analysis_scope,
  reversible         = reversible,
  total_clones       = nrow(trees),
  converged_clones   = nrow(trees_converged),
  unconverged_clones = length(unconverged_idx),
  unconverged_clone_ids = if (length(unconverged_idx) > 0) paste(unconverged_clone_ids, collapse = ";") else "",
  final_iteration    = max_iter,
  mcmc_length        = DEFAULT_MCMC_LENGTH,
  ess_cutoff         = ess_cutoff,
  ignore_params      = paste(ignore_params, collapse = ";"),
  timestamp          = as.character(Sys.time()),
  stringsAsFactors   = FALSE
)

# Save individual convergence summary
individual_summary_file <- file.path(summary_dir, paste0(config_name, "_", template_id, "_convergence_summary.csv"))
write.csv(convergence_summary, individual_summary_file, row.names = FALSE)
cat("Individual convergence summary saved to:", individual_summary_file, "\n")

# Update combined convergence summary
combined_summary_file <- file.path(summary_dir, paste0(model_type, "_convergence_summary.csv"))

if (file.exists(combined_summary_file)) {
  # Read existing summary and update
  existing_summary <- read.csv(combined_summary_file, stringsAsFactors = FALSE)

  # Remove any existing entry for this config/template combination
  existing_summary <- existing_summary[!(existing_summary$config_name == config_name &
                                           existing_summary$template_id == template_id), ]

  # Add new entry
  updated_summary <- rbind(existing_summary, convergence_summary)
  write.csv(updated_summary, combined_summary_file, row.names = FALSE)
  cat("Updated combined convergence summary:", combined_summary_file, "\n")
} else {
  # Create new combined summary
  write.csv(convergence_summary, combined_summary_file, row.names = FALSE)
  cat("Created new combined convergence summary:", combined_summary_file, "\n")
}

# ------------------ GC clock rates estimation ------------------ #
# Post-processing: extract individual clone clock rates for StrictClock phases
if (model_type == "gc_strict_clock") {
  cat("=== Extracting GC clock rates ===\n")
  SC_trees <- readRDS(converged_trees_file)
  gc_clock_rates <- sapply(SC_trees$parameters, function(x) {
    x$mean[x$item == "geneticClockRate"]
  })

  if (sum(is.na(gc_clock_rates)) > 0) {
    cat("Warning: There are NA clock rates calculated for the GC B cells.")
    num_nas <- sum(is.na(gc_clock_rates))
    cat("Number of NA clock rates:", num_nas)
  }

  mean_gc_clock_rate <- mean(gc_clock_rates, na.rm = TRUE)
  individual_mean_rates_file <- file.path(clock_rates_dir, paste0(config_name, "_individual_mean_gc_clock_rates.csv"))

  new_mean_entry <- data.frame(
    config_name     = config_name,
    mean_clock_rate = mean_gc_clock_rate,
    num_clones      = length(gc_clock_rates[!is.na(gc_clock_rates)]),
    timestamp       = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.csv(new_mean_entry, individual_mean_rates_file)

  cat("Individual mean clone clock rates saved to:", individual_mean_rates_file, "\n")

  if (length(gc_clock_rates) > 0) {
    cat("Mean:", round(mean_gc_clock_rate, 8), "\n")
    cat("Median:", round(median(gc_clock_rates, na.rm = TRUE), 8), "\n")
    cat("Range:", round(min(gc_clock_rates, na.rm = TRUE), 8), "-",
        round(max(gc_clock_rates, na.rm = TRUE), 8), "\n")
  }
}

# ------------------ Final logging and cleanup ------------------ #
cat("=== Analysis Summary ===\n")
logKVs("Final Results", c(
  Total_Clones       = nrow(trees),
  Converged_Clones   = nrow(trees_converged),
  Convergence_Rate   = paste0(round(100 * nrow(trees_converged) / nrow(trees), 1), "%"),
  All_Trees_File     = all_trees_file,
  Converged_Trees_File = converged_trees_file,
  Beast_Output_Dir   = beast_raw_output_dir,
  Summary_File       = individual_summary_file
))

cat("=== BEAST Dowser Job", task_id, "Completed at", as.character(Sys.time()), "===\n")