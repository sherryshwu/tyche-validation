#!/bin/bash
# Get parameters
SIMULATION_NAME="${1:-tltt_08_20}"
ANALYSIS_TYPE="${2:-main_analysis}"
REV_SUFFIX="${3:-irrev}"

# Setup paths
PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
RESULTS_BASE_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/results/${ANALYSIS_TYPE}/${REV_SUFFIX}"
TREE_ANALYSIS_DIR="${RESULTS_BASE_DIR}/tree_analysis"
RAW_DATA_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/data/raw"
LOG_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/logs/${ANALYSIS_TYPE}/tree_analysis"

# Setup logging
LOG_FILE="${LOG_DIR}/tree_analysis_${ANALYSIS_TYPE}_${REV_SUFFIX}_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG_FILE")"

# Redirect all output to both console and log file
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "=== Tree Analysis Pipeline Started: $(date) ==="
echo "Simulation: $SIMULATION_NAME"
echo "Analysis: $ANALYSIS_TYPE ($REV_SUFFIX)"
echo "Log file: $LOG_FILE"

# Create directories
mkdir -p "$TREE_ANALYSIS_DIR" "$LOG_DIR"

echo "Results directory: $RESULTS_BASE_DIR"
echo "Tree analysis: $TREE_ANALYSIS_DIR"

# Setup environment
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate r_phylo
cd "$PROJECT_ROOT"

# =============================================================================
# Step 1: Find all beast outputs and create job lists
echo ""
echo "=== Step 1: Finding BEAST Outputs ==="

# Setup configs directory
CONFIGS_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/configs"
mkdir -p "$CONFIGS_DIR"

# Find all model types that have been run
MODEL_TYPES=()
for model_dir in "$RESULTS_BASE_DIR"/*; do
    if [[ -d "$model_dir" ]]; then
        model_type=$(basename "$model_dir")
        MODEL_TYPES+=("$model_type")
        echo "Found model type: $model_type"
    fi
done

if [[ ${#MODEL_TYPES[@]} -eq 0 ]]; then
    echo "ERROR: No model directories found in $RESULTS_BASE_DIR"
    exit 1
fi

echo "Total model types found: ${#MODEL_TYPES[@]}"

# Create CSV job list in configs directory
JOB_LIST_FILE="$CONFIGS_DIR/tree_analysis_jobs_${ANALYSIS_TYPE}_${REV_SUFFIX}.csv"

# Write CSV header
echo "job_id,config_name,template_name,beast_tree_files,true_tree_file,model_type,analysis_type,rev_suffix" > "$JOB_LIST_FILE"

TOTAL_JOBS=0

for model_type in "${MODEL_TYPES[@]}"; do
    echo "Processing model type: $model_type"
    
    # Find all beast_raw_output directories
    beast_output_dirs=($(find "$RESULTS_BASE_DIR/$model_type" -name "beast_raw_output" -type d))
    
    if [[ ${#beast_output_dirs[@]} -eq 0 ]]; then
        echo "  Warning: No beast_raw_output found for $model_type"
        continue
    fi
    
    for beast_output_dir in "${beast_output_dirs[@]}"; do
        echo "  Scanning: $beast_output_dir"
        
        # Add timeout and error handling
        config_paths=($(timeout 60 find "$beast_output_dir" -maxdepth 1 -type d -not -path "$beast_output_dir" 2>/dev/null))
        
        if [[ $? -ne 0 ]]; then
            echo "    ERROR: Timeout or error scanning $beast_output_dir"
            continue
        fi

        echo "    Found ${#config_paths[@]} config directories"

        for config_path in "$beast_output_dir"/*; do
            if [[ -d "$config_path" ]]; then
                config_name=$(basename "$config_path")
                
                # Find true tree file
                true_tree_file="$RAW_DATA_DIR/$config_name/all_simplified_time_trees.nex"
                if [[ ! -f "$true_tree_file" ]]; then
                    echo "    Warning: True tree not found for $config_name"
                    continue
                fi
                
                # Find BEAST tree files
                beast_tree_files=($(find "$config_path" -name "*_tree_with_trait*.tree" | sort))
                if [[ $? -ne 0 ]]; then
                    echo "      ERROR: Timeout finding tree files in $config_path"
                    continue
                fi
                if [[ ${#beast_tree_files[@]} -eq 0 ]]; then
                    echo "    Warning: No BEAST trees found for $config_name"
                    continue
                fi
                
                echo "      Found ${#beast_tree_files[@]} tree files"
                echo "    Grouping templates..."

                # Group by template
                templates=($(printf '%s\n' "${beast_tree_files[@]}" | \
                            sed 's|.*/||; s/_[0-9]\+_tree_with_trait.*\.tree$//' | \
                            sort -u))
                
                echo "    Found ${#templates[@]} templates: ${templates[*]}"
                
                for template_name in "${templates[@]}"; do
                    echo "      Processing template: $template_name"
                    # Get all files for this template
                    beast_tree_file=($(printf '%s\n' "${beast_tree_files[@]}" | \
                                    grep "/${template_name}_[0-9]\+_tree_with_trait"))
                    
                    # Convert array to semicolon-separated string
                    beast_tree_files=$(IFS=';'; echo "${beast_tree_file[*]}")
                    
                    # Increment job counter
                    ((TOTAL_JOBS++))
                    
                    # Add to CSV with additional metadata
                    echo "$TOTAL_JOBS,\"$config_name\",\"$template_name\",\"$beast_tree_files\",\"$true_tree_file\",\"$model_type\",\"$ANALYSIS_TYPE\",\"$REV_SUFFIX\"" >> "$JOB_LIST_FILE"
                done
                
                echo "    Added jobs for $config_name: ${#templates[@]} templates"
            fi
        done
    done
done

echo "Total jobs created: $TOTAL_JOBS"
echo "Job list saved to: $JOB_LIST_FILE"

# =============================================================================
# Step 1.5: Generate convergence summaries if missing
echo ""
echo "=== Step 1.5: Checking Convergence Summaries ==="

# Check if convergence summaries exist
TYCHE_CONV_FILE="$RESULTS_BASE_DIR/tyche_models/summary/tyche_models_convergence_summary.csv"
COMP_CONV_FILE="$RESULTS_BASE_DIR/competing_models/summary/competing_models_convergence_summary.csv"

if [[ ! -f "$TYCHE_CONV_FILE" ]] || [[ ! -f "$COMP_CONV_FILE" ]]; then
    echo "Convergence summaries missing - generating them..."
    
    Rscript scripts/analysis/summarize_convergence.R \
        "$SIMULATION_NAME" "$REV_SUFFIX" "$ANALYSIS_TYPE"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Convergence summaries generated successfully"
    else
        echo "ERROR: Failed to generate convergence summaries"
        exit 1
    fi
else
    echo "✓ Convergence summaries already exist"
fi

# =============================================================================
# Step 2: Submit Slurm jobs
echo ""
echo "=== Step 2: Submitting Analysis Jobs ==="

if [[ $TOTAL_JOBS -eq 0 ]]; then
    echo "ERROR: No analysis jobs to submit"
    exit 1
fi

# Verify job list exists
if [[ ! -f "$JOB_LIST_FILE" ]]; then
    echo "ERROR: Job list file not found: $JOB_LIST_FILE"
    exit 1
fi

echo "Using job list: $JOB_LIST_FILE"

ANALYSIS_LOG_DIR="${LOG_DIR}/tree_analysis_jobs"
mkdir -p "$ANALYSIS_LOG_DIR"

JOB_ID=$(sbatch --parsable \
    --array=1-${TOTAL_JOBS}%10 \
    --cpus-per-task=4 \
    --mem-per-cpu=6gb \
    --time=2:00:00 \
    --job-name=tree-analysis-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="${ANALYSIS_LOG_DIR}/tree_analysis_%A_%a.out" \
    --error="${ANALYSIS_LOG_DIR}/tree_analysis_%A_%a.err" \
    --account=hoehnlab-share \
    --wrap="
        echo \"Job \$SLURM_ARRAY_TASK_ID started: \$(date)\"
        source /optnfs/common/miniconda3/etc/profile.d/conda.sh
        conda activate r_phylo
        cd \"$PROJECT_ROOT\"
        
        Rscript scripts/analysis/calculate_tree_metrics.R \
            \"$JOB_LIST_FILE\" \
            \"\$SLURM_ARRAY_TASK_ID\" \
            \"$TREE_ANALYSIS_DIR\" \
            \"$SIMULATION_NAME\" \
            \"$ANALYSIS_TYPE\" \
            \"$REV_SUFFIX\"
        
        # Log job completion
        echo \"Job \$SLURM_ARRAY_TASK_ID completed: \$(date)\"
    ")

echo "Jobs submitted with ID: $JOB_ID"
echo "Job list: $JOB_LIST_FILE"
echo "Waiting for analysis jobs to complete before proceeding..."

# =============================================================================
# Step 3: Job status check
echo ""
echo "=== Step 3: Job Status Check ==="
sleep 5  # Brief pause to let jobs start

job_status=$(squeue -j "$JOB_ID" 2>/dev/null | wc -l)
if [[ $job_status -gt 1 ]]; then
    running_jobs=$((job_status - 1))
    echo "✓ $running_jobs jobs are running (Job ID: $JOB_ID)"
else
    echo "No jobs found in queue - they may have finished very quickly or failed to start"
fi

# =============================================================================
# Step 4: Create combined summary (depends on analysis jobs)
echo ""
echo "=== Step 4: Creating Combined Summary ===" 
SUMMARY_LOG="${LOG_DIR}/summary_job.log"
SUMMARY_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${JOB_ID} \
    --cpus-per-task=1 \
    --mem-per-cpu=2gb \
    --time=10:00 \
    --job-name=tree-analysis-summary-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="${SUMMARY_LOG}" \
    --error="${SUMMARY_LOG}" \
    --account=hoehnlab-share \
    --wrap="
    echo \"=== Summary job started at \$(date) ===\" 
    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
    conda activate r_phylo
    cd \"$PROJECT_ROOT\"
    python3 scripts/analysis/consolidate_tree_metrics.py \"$TREE_ANALYSIS_DIR\"
    echo \"=== Summary job completed at \$(date) ===\" 
    ")
echo "Creating combined summary job submitted with ID: $SUMMARY_JOB_ID"

# =============================================================================
# Step 5: Submit Tree Plotting Jobs
echo ""
echo "=== Step 5: Submitting Tree Plotting Jobs ==="
TREE_PLOTS_LOG_DIR="${LOG_DIR}/tree_plotting_jobs"
mkdir -p "$TREE_PLOTS_LOG_DIR"

PLOT_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${SUMMARY_JOB_ID} \
    --array=1-${TOTAL_JOBS}%5 \
    --cpus-per-task=2 \
    --mem-per-cpu=8gb \
    --time=1:00:00 \
    --job-name=tree-plots-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="${TREE_PLOTS_LOG_DIR}/tree_plots_%A_%a.out" \
    --error="${TREE_PLOTS_LOG_DIR}/tree_plots_%A_%a.err" \
    --account=hoehnlab-share \
    --wrap="
        echo \"Tree plotting job \$SLURM_ARRAY_TASK_ID started: \$(date)\"
        source /optnfs/common/miniconda3/etc/profile.d/conda.sh
        conda activate r_phylo
        cd \"$PROJECT_ROOT\"
        Rscript scripts/plotting/plot_trees.R \
            \"$JOB_LIST_FILE\" \
            \"\$SLURM_ARRAY_TASK_ID\" \
            \"$TREE_ANALYSIS_DIR\" \
            \"$ANALYSIS_TYPE\"
        
        echo \"Tree plotting job \$SLURM_ARRAY_TASK_ID completed: \$(date)\"
    ")

echo "Tree plotting jobs submitted with ID: $PLOT_JOB_ID"
echo "Individual job logs in: $TREE_PLOTS_LOG_DIR"

echo ""
echo "=== Tree Analysis Pipeline Completed: $(date) ==="
echo "Main log file: $LOG_FILE"
echo "Analysis job logs: $ANALYSIS_LOG_DIR/"
echo "Tree plotting logs: $TREE_PLOTS_LOG_DIR/"
echo "Results directory: $TREE_ANALYSIS_DIR"