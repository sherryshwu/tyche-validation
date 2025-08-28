#!/bin/bash
# =============================================================================
# TREE ANALYSIS PIPELINE - STREAMLINED FOR PUBLICATION
# =============================================================================

set -e

# Get parameters
SIMULATION_NAME="${1:-tltt_08_20}"
ANALYSIS_TYPE="${2:-main_analysis}"
REV_SUFFIX="${3:-rev}"

# Setup paths
PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
RESULTS_BASE_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/results/${ANALYSIS_TYPE}/${REV_SUFFIX}"
TREE_ANALYSIS_DIR="${RESULTS_BASE_DIR}/tree_analysis"
RAW_DATA_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/data/raw"
LOG_DIR="${PROJECT_ROOT}/${SIMULATION_NAME}/logs/${ANALYSIS_TYPE}"

# Setup logging
LOG_FILE="${PROJECT_ROOT}/${SIMULATION_NAME}/logs/${ANALYSIS_TYPE}/tree_analysis_${ANALYSIS_TYPE}_${REV_SUFFIX}_$(date +%Y%m%d_%H%M%S).log"
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

#-----------------------------------------------------
# Step 1: Find all beast outputs and create job lists
#-----------------------------------------------------
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
echo "job_id,config_name,template_name,beast_files,true_tree_file,model_type,analysis_type,rev_suffix" > "$JOB_LIST_FILE"

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
                beast_files=($(find "$config_path" -name "*_tree_with_trait*.tree" | sort))
                if [[ $? -ne 0 ]]; then
                    echo "      ERROR: Timeout finding tree files in $config_path"
                    continue
                fi
                if [[ ${#beast_files[@]} -eq 0 ]]; then
                    echo "    Warning: No BEAST trees found for $config_name"
                    continue
                fi
                
                echo "      Found ${#beast_files[@]} tree files"

                # Group by template
                templates=($(printf '%s\n' "${beast_files[@]}" | \
                            sed 's|.*/||; s/_[0-9]\+_tree_with_trait.*\.tree$//' | \
                            sort -u))
                
                for template_name in "${templates[@]}"; do
                    # Get all files for this template
                    beast_tree_file=($(printf '%s\n' "${beast_files[@]}" | \
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

JOB_ID=$(sbatch --parsable \
    --array=1-${TOTAL_JOBS}%10 \
    --cpus-per-task=4 \
    --mem-per-cpu=6gb \
    --time=2:00:00 \
    --job-name=tree-analysis-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out" \
    --error="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err" \
    --account=hoehnlab-share \
    --wrap="
        source /optnfs/common/miniconda3/etc/profile.d/conda.sh
        conda activate r_phylo
        cd \"$PROJECT_ROOT\"
        
        Rscript scripts/tree_analysis.R \
            \"$JOB_LIST_FILE\" \
            \"\$SLURM_ARRAY_TASK_ID\" \
            \"$TREE_ANALYSIS_DIR\"
    ")

echo "Jobs submitted with ID: $JOB_ID"
echo "Job list: $JOB_LIST_FILE"

# Save job ID for monitoring
echo "$JOB_ID" > "$LOG_DIR/tree_analysis_job_id.txt"
echo "Job ID saved to: $LOG_DIR/tree_analysis_job_id.txt"

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
# Step 4: Create combined summary
echo ""
echo "=== Step 4: Creating Combined Summary ===" 
python3 scripts/create_combined_summary.py "$TREE_ANALYSIS_DIR"

# =============================================================================
# Step 5: Creating publication plots
echo ""
echo "=== Step 5: Creating Publication Plots ==="
selected_models="EO_Fixed,EO_Est,IS_Est,MS_Fixed,MS_Est,SC_AR,UCLD_AR"
Rscript scripts/publication_plots.R "$TREE_ANALYSIS_DIR" "$SIMULATION_NAME" "$ANALYSIS_TYPE" "$REV_SUFFIX" "$selected_models"

# =============================================================================
# Step 6: Submit Tree Plotting Jobs
echo ""
echo "=== Step 6: Submitting Tree Plotting Jobs ==="

PLOT_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${JOB_ID} \
    --array=1-${TOTAL_JOBS}%5 \
    --cpus-per-task=2 \
    --mem-per-cpu=8gb \
    --time=1:00:00 \
    --job-name=tree-plots-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out" \
    --error="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err" \
    --account=hoehnlab-share \
    --wrap="
        source /optnfs/common/miniconda3/etc/profile.d/conda.sh
        conda activate r_phylo
        cd \"$PROJECT_ROOT\"
        
        # Tree plotting
        Rscript scripts/tree_plotting.R \
            \"$JOB_LIST_FILE\" \
            \"\$SLURM_ARRAY_TASK_ID\" \
            \"$TREE_ANALYSIS_DIR\"
    ")
# sbatch --array=1-${TOTAL_JOBS}%5 \
#     --cpus-per-task=2 \
#     --mem-per-cpu=8gb \
#     --time=1:00:00 \
#     --job-name=tree-plots-${ANALYSIS_TYPE}-${REV_SUFFIX} \
#     --output="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out" \
#     --error="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err" \
#     --account=hoehnlab-share \
#     --wrap="
#         source /optnfs/common/miniconda3/etc/profile.d/conda.sh
#         conda activate r_phylo
#         cd \"$PROJECT_ROOT\"
        
#         # Tree plotting
#         Rscript scripts/tree_plotting.R \
#             \"$JOB_LIST_FILE\" \
#             \"\$SLURM_ARRAY_TASK_ID\" \
#             \"$TREE_ANALYSIS_DIR\"
#     "

echo "Tree plotting jobs submitted with ID: $PLOT_JOB_ID"
echo "Plot job ID saved to: $LOG_DIR/tree_plotting_job_id.txt"
echo "$PLOT_JOB_ID" > "$LOG_DIR/tree_plotting_job_id.txt"

# =============================================================================
# Step 7: Submit Publication Aligned Trees Jobs
echo ""
echo "=== Step 7: Submitting Publication Aligned Trees Jobs ==="

PUB_TREES_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${PLOT_JOB_ID} \
    --array=1-${TOTAL_JOBS}%10 \
    --cpus-per-task=2 \
    --mem-per-cpu=4gb \
    --time=30:00 \
    --job-name=pub-trees-${ANALYSIS_TYPE}-${REV_SUFFIX} \
    --output="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out" \
    --error="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err" \
    --account=hoehnlab-share \
    --wrap="
        source /optnfs/common/miniconda3/etc/profile.d/conda.sh
        conda activate r_phylo
        cd \"$PROJECT_ROOT\"
        
        # Publication aligned trees (default to clone 3)
        Rscript scripts/plot_aligned_trees.R \
            \"$JOB_LIST_FILE\" \
            \"\$SLURM_ARRAY_TASK_ID\" \
            \"$TREE_ANALYSIS_DIR\" \
            \"3\"
    ")

# PUB_TREES_JOB_ID=$(sbatch --parsable \
#     --array=1-${TOTAL_JOBS}%10 \
#     --cpus-per-task=2 \
#     --mem-per-cpu=4gb \
#     --time=30:00 \
#     --job-name=pub-trees-${ANALYSIS_TYPE}-${REV_SUFFIX} \
#     --output="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out" \
#     --error="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err" \
#     --account=hoehnlab-share \
#     --wrap="
#         source /optnfs/common/miniconda3/etc/profile.d/conda.sh
#         conda activate r_phylo
#         cd \"$PROJECT_ROOT\"
        
#         # Publication aligned trees (default to clone 3)
#         Rscript scripts/plot_scaled_trees.R \
#             \"$JOB_LIST_FILE\" \
#             \"\$SLURM_ARRAY_TASK_ID\" \
#             \"$TREE_ANALYSIS_DIR\" \
#             \"3\"
#     ")

echo "Publication aligned trees jobs submitted with ID: $PUB_TREES_JOB_ID"
echo "Pub trees job ID saved to: $LOG_DIR/pub_trees_job_id.txt"
echo "$PUB_TREES_JOB_ID" > "$LOG_DIR/pub_trees_job_id.txt"

# Update the final echo statement
echo ""
echo "=== Tree Analysis Pipeline Completed: $(date) ==="
echo "Results directory: $TREE_ANALYSIS_DIR"
echo "Check publication plots in: $TREE_ANALYSIS_DIR/plots/publication_plots"
echo "Check aligned trees in: $TREE_ANALYSIS_DIR/plots/publication_plots/aligned_trees"