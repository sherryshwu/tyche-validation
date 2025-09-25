#!/bin/bash
#SBATCH --job-name=run-beast-dowser
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=144:00:00
#SBATCH --output=/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out
#SBATCH --error=/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err
#SBATCH --account=hoehnlab
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

set -e

# Get parameters from command line
simulation_run=${1:-"tltt_08_20"}
analysis_scope=${2}                  # Required: main_analysis, sub_analysis, differentiation_analysis
model_type=${3}                      # Required: gc_strict_clock, tyche_models, competing_models
reversible=${4:-false}               # Whether non-GC to GC transitions are allowed
nproc=${5:-20}                       # Number of processors for BEAST

# Define run parameters
PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
SCRATCH_ROOT="/dartfs-hpc/scratch/f0070d5"
BEAST_DIR="/dartfs/rc/lab/H/HoehnK/Sherry/beast/bin"
DEFAULT_SEED=12345

# Determine reversibility suffix for directory naming
if [ "$reversible" = "true" ]; then
    rev_suffix="rev"
else
    rev_suffix="irrev"
fi

# Define directory structure
BASE_DATA_DIR="${PROJECT_ROOT}/${simulation_run}/data"
RAW_DATA_DIR="${BASE_DATA_DIR}/raw"
PROCESSED_DATA_DIR="${BASE_DATA_DIR}/processed"
ANALYSIS_PROCESSED_DATA_DIR="${PROCESSED_DATA_DIR}/${analysis_scope}"
DOWSER_PROCESSED_TREES_DIR="${ANALYSIS_PROCESSED_DATA_DIR}/dowser_processed_trees/${rev_suffix}"
CONFIGS_DIR="${PROJECT_ROOT}/${simulation_run}/configs"

RESULTS_BASE_DIR="${PROJECT_ROOT}/${simulation_run}/results/${analysis_scope}/${rev_suffix}/${model_type}"

BEAST_RAW_OUTPUT_DIR="${RESULTS_BASE_DIR}/beast_raw_output"
SUMMARY_DIR="${RESULTS_BASE_DIR}/summary"
CORRELATION_DIR="${RESULTS_BASE_DIR}/correlation_analysis"
CLOCK_RATES_DIR="${RESULTS_BASE_DIR}/clock_rates"

# Scratch directory for temporary BEAST output
SCRATCH_BEAST_DIR="${SCRATCH_ROOT}/${simulation_run}/${analysis_scope}/${rev_suffix}/${model_type}"

# Centralized logs directory
LOG_DIR="${PROJECT_ROOT}/${simulation_run}/logs/${analysis_scope}/beast_runs"

# Create necessary directories
mkdir -p "$RAW_DATA_DIR" "$PROCESSED_DATA_DIR" "$ANALYSIS_PROCESSED_DATA_DIR" "$DOWSER_PROCESSED_TREES_DIR" "$CONFIGS_DIR"
mkdir -p "$BEAST_RAW_OUTPUT_DIR" "$SUMMARY_DIR" "$SCRATCH_BEAST_DIR" "$LOG_DIR"

if [[ "$model_type" == "gc_strict_clock" ]]; then
    mkdir -p "$CORRELATION_DIR" "$CLOCK_RATES_DIR"
fi

# Load conda environment
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda --no-plugins init bash
conda activate r_phylo

# Change to project directory
cd "$PROJECT_ROOT"

# Get current task ID from SLURM array
TASK_ID=${SLURM_ARRAY_TASK_ID}

# Determine job combinations file based on analysis scope and model type
determine_JOB_COMBINATIONS_FILE() {
    local base_name="beast_job_combinations"
    local suffix="${model_type}.csv"
    
    if [[ "$analysis_scope" == "main_analysis" ]]; then
        echo "${CONFIGS_DIR}/${base_name}_${suffix}"
    else
        echo "${CONFIGS_DIR}/${base_name}_${analysis_scope}_${suffix}"
    fi
}

JOB_COMBINATIONS_FILE=$(determine_JOB_COMBINATIONS_FILE)

# Define step names for logging
case "$model_type" in
    gc_strict_clock)
        step_description="Strict Clock on Germinal Center B cells"
        ;;
    tyche_models)
        step_description="TyCHE models (5 variants) on all B cell populations"
        ;;
    competing_models)
        step_description="Competing models (Strict + UCLD) on all B cell populations"
        ;;
esac

# Extract config name from combinations file for this task
config_name=$(sed -n "${TASK_ID}p" <(tail -n +2 "$JOB_COMBINATIONS_FILE") | cut -d',' -f1)
template_id=$(sed -n "${TASK_ID}p" <(tail -n +2 "$JOB_COMBINATIONS_FILE") | cut -d',' -f2)

# Create analysis-specific log file
if [[ "$analysis_scope" == "sub_analysis" ]]; then
    time_subset=$(sed -n "${TASK_ID}p" <(tail -n +2 "$JOB_COMBINATIONS_FILE") | cut -d',' -f5)
    if [[ -n "$time_subset" ]]; then
        LOG_FILE="${LOG_DIR}/beast_dowser_${config_name}_${template_id}_time${time_subset}_${rev_suffix}.log"
    else
        LOG_FILE="${LOG_DIR}/beast_dowser_${config_name}_${template_id}_${rev_suffix}.log"
    fi
else
    LOG_FILE="${LOG_DIR}/beast_dowser_${config_name}_${template_id}_${rev_suffix}.log"
fi

echo "=== Starting BEAST dowser analysis - $step_description at $(date) ===" | tee "$LOG_FILE"
echo "Task ID: $TASK_ID" | tee -a "$LOG_FILE"
echo "Simulation run: $simulation_run" | tee -a "$LOG_FILE"
echo "Analysis scope: $analysis_scope" | tee -a "$LOG_FILE"
echo "Model Type: $model_type ($step_description)" | tee -a "$LOG_FILE"
echo "Reversible Transitions: $reversible (directory suffix: $rev_suffix)" | tee -a "$LOG_FILE"
echo "Number of processes: $nproc" | tee -a "$LOG_FILE"
echo "Job combinations file: $JOB_COMBINATIONS_FILE" | tee -a "$LOG_FILE"
echo "Permanent Beast output directory: $BEAST_RAW_OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Scratch Beast output directory: $SCRATCH_BEAST_DIR" | tee -a "$LOG_FILE"
echo "Dowser Processed Trees: $DOWSER_PROCESSED_TREES_DIR" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

# Verify job combinations file exists
if [[ ! -f "$JOB_COMBINATIONS_FILE" ]]; then
    echo "ERROR: Job combinations file not found: $JOB_COMBINATIONS_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

# Verify task ID is valid for this phase
total_jobs=$(tail -n +2 "$JOB_COMBINATIONS_FILE" | wc -l)
if [[ $TASK_ID -gt $total_jobs ]]; then
    echo "ERROR: Task ID $TASK_ID exceeds total jobs ($total_jobs) for $step_description" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Running BEAST dowser analysis - $step_description..." | tee -a "$LOG_FILE"

# Run single BEAST job with all parameters
Rscript scripts/phylogenetics/run_beast_dowser.R \
"${TASK_ID}" \
"$JOB_COMBINATIONS_FILE" \
"$simulation_run" \
"$analysis_scope" \
"$model_type" \
"$reversible" \
"$BEAST_DIR" \
"$nproc" \
"$DEFAULT_SEED" 2>&1 | tee -a "$LOG_FILE"

echo "=== BEAST dowser analysis - $step_description completed at $(date) ===" | tee -a "$LOG_FILE"

# Copy results from scratch to permanent storage
echo "=== Starting results copy at $(date) ===" | tee -a "$LOG_FILE"

# Updated paths with phase subdirectory
if [[ "$analysis_scope" == "sub_analysis" ]]; then
    # Extract time subset from the job combinations file (5th column)
    time_subset=$(sed -n "${TASK_ID}p" <(tail -n +2 "$JOB_COMBINATIONS_FILE") | cut -d',' -f5)
    
    if [[ -z "$time_subset" ]]; then
        echo "ERROR: No time subset found for sub_analysis task $TASK_ID" | tee -a "$LOG_FILE"
        echo "Expected time subset in column 5 of job combinations file" | tee -a "$LOG_FILE"
        exit 1
    fi
    
    echo "Time subset: $time_subset" | tee -a "$LOG_FILE"
    
    # Update paths to include time subset
    scratch_output_dir="${SCRATCH_BEAST_DIR}/time_${time_subset}/${config_name}"
    perm_output_dir="${BEAST_RAW_OUTPUT_DIR}/time_${time_subset}/${config_name}"
else
    # Standard paths for main_analysis and differentiation_analysis
    scratch_output_dir="${SCRATCH_BEAST_DIR}/${config_name}"
    perm_output_dir="${BEAST_RAW_OUTPUT_DIR}/${config_name}"
fi

echo "Config name: $config_name" | tee -a "$LOG_FILE"
echo "Scratch directory: $scratch_output_dir" | tee -a "$LOG_FILE"
echo "Permanent directory: $perm_output_dir" | tee -a "$LOG_FILE"

# Copy results if scratch directory exists
if [[ -d "$scratch_output_dir" ]]; then
    # Create permanent directory
    mkdir -p "$perm_output_dir"
    
    # Copy all files from scratch directory
    if cp -r "$scratch_output_dir"/* "$perm_output_dir"/; then
        echo "Results successfully copied" | tee -a "$LOG_FILE"

    else
        echo "ERROR: Failed to copy results" | tee -a "$LOG_FILE"
        exit 1
    fi
else
    echo "WARNING: Scratch directory not found: $scratch_output_dir" | tee -a "$LOG_FILE"
fi

echo "=== Results copy completed at $(date) ===" | tee -a "$LOG_FILE"
echo "=== Analysis pipeline completed successfully at $(date) ===" | tee -a "$LOG_FILE"