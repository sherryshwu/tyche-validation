#!/bin/bash
# Multi-analysis publication plotting script

# Get parameters
ANALYSIS_TYPE="${1:-main_analysis}"
BASE_PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"

echo "=== Creating Publication Plots: $(date) ==="
echo "Analysis type: $ANALYSIS_TYPE"

# Setup combined output directory
COMBINED_PLOTS_DIR="${BASE_PROJECT_ROOT}/figures/${ANALYSIS_TYPE}"
mkdir -p "$COMBINED_PLOTS_DIR"

# Setup logging
LOG_FILE="${COMBINED_PLOTS_DIR}/publication_plots_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE")
exec 2>&1

# =============================================================================
# Step 1: Main Publication Plots (Cross-simulation comparison)
echo ""
echo "=== Step 1: Creating Main Publication Plots ==="
echo "Target: Two simulations x Two 1:1 configs = Four panels"

selected_models="EO_Fixed,EO_Est,IS_Est,SC_AR,UCLD_AR"
selected_configs="config_ratio_1to1_sel,config_ratio_1to1_neu"
selected_simulations="tltt_08_20,gc_reentry_hunter"
selected_rev_suffixes="irrev"

# Build summary file paths
SUMMARY_FILES=""
for sim in tltt_08_20 gc_reentry_hunter; do
    for rev in irrev; do
        summary_file="${BASE_PROJECT_ROOT}/${sim}/results/${ANALYSIS_TYPE}/${rev}/tree_analysis/all_results_summary.csv"
        if [[ -f "$summary_file" ]]; then
            if [[ -n "$SUMMARY_FILES" ]]; then
                SUMMARY_FILES="${SUMMARY_FILES},${summary_file}"
            else
                SUMMARY_FILES="$summary_file"
            fi
        else
            echo "Warning: Summary file not found: $summary_file"
        fi
    done
done

echo "Summary files: $SUMMARY_FILES"

MAIN_PLOTS_JOB_ID=$(sbatch --parsable \
    --cpus-per-task=2 \
    --mem-per-cpu=6gb \
    --time=45:00 \
    --job-name=main-publication-plots \
    --output="${LOG_FILE}" \
    --error="${LOG_FILE}" \
    --account=hoehnlab-share \
    --wrap="
    echo \"Main publication plots started: \$(date)\"
    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
    conda activate r_phylo
    cd \"$BASE_PROJECT_ROOT\"
    
    Rscript scripts/plotting/publication_plots_main.R \
        \"$SUMMARY_FILES\" \
        \"$selected_simulations\" \
        \"$ANALYSIS_TYPE\" \
        \"$selected_rev_suffixes\" \
        \"$selected_models\" \
        \"$selected_configs\" \
        \"main\" \
        \"$COMBINED_PLOTS_DIR\"
    
    echo \"Main publication plots completed: \$(date)\"
    ")

# =============================================================================
# Step 2: Supplementary Plots (Single simulation, multiple configs)
echo ""
echo "=== Step 2: Creating Supplementary Publication Plots ==="
echo "Target: One simulation x Two configs (1:3, sel and neu) = Two panels"

selected_configs_supp="config_ratio_1to3_sel,config_ratio_1to3_neu,config_ratio_1to1_sel,config_ratio_1to1_neu"
selected_simulations_supp="tltt_08_20"

# Build summary file for supplementary
SUMMARY_FILE_SUPP="${BASE_PROJECT_ROOT}/tltt_08_20/results/${ANALYSIS_TYPE}/irrev/tree_analysis/all_results_summary.csv"

SUP_PLOTS_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${MAIN_PLOTS_JOB_ID} \
    --cpus-per-task=2 \
    --mem-per-cpu=6gb \
    --time=45:00 \
    --job-name=supp-publication-plots \
    --output="${LOG_FILE}" \
    --error="${LOG_FILE}" \
    --account=hoehnlab-share \
    --wrap="
    echo \"Supplementary plots started: \$(date)\"
    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
    conda activate r_phylo
    cd \"$BASE_PROJECT_ROOT\"
    
    Rscript scripts/plotting/publication_plots_main.R \
        \"$SUMMARY_FILE_SUPP\" \
        \"$selected_simulations_supp\" \
        \"$ANALYSIS_TYPE\" \
        \"irrev\" \
        \"$selected_models\" \
        \"$selected_configs_supp\" \
        \"supplementary\" \
        \"$COMBINED_PLOTS_DIR\"
    
    echo \"Supplementary plots completed: \$(date)\"
    ")

echo ""
echo "=== Publication Plots Pipeline Submitted ==="
echo "Main plots job: $MAIN_PLOTS_JOB_ID"
echo "Supplementary plots job: $SUP_PLOTS_JOB_ID"
echo "Output directory: $COMBINED_PLOTS_DIR"
echo "Log file: $LOG_FILE"