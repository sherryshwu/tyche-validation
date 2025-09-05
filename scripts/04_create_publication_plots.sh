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

# Common settings
selected_models="EO_Fixed,EO_Est,IS_Est,SC_AR,UCLD_AR"
selected_rev_suffixes="irrev"

# =============================================================================
# Step 1: Main Publication Plots (Cross-simulation comparison)
echo ""
echo "=== Step 1: Creating Main Publication Plots ==="
echo "Target: Two simulations x Two 1:1 configs = Four panels"

selected_configs_main="config_ratio_1to1_sel,config_ratio_1to1_neu"
selected_simulations="tltt_08_20,gc_reentry_hunter"

# Build summary file paths for both simulations
SUMMARY_FILES_BOTH=""
for sim in tltt_08_20 gc_reentry_hunter; do
    summary_file="${BASE_PROJECT_ROOT}/${sim}/results/${ANALYSIS_TYPE}/irrev/tree_analysis/all_results_summary.csv"
    if [[ -f "$summary_file" ]]; then
        if [[ -n "$SUMMARY_FILES_BOTH" ]]; then
            SUMMARY_FILES_BOTH="${SUMMARY_FILES_BOTH},${summary_file}"
        else
            SUMMARY_FILES_BOTH="$summary_file"
        fi
    fi
done

echo "Summary files for main: $SUMMARY_FILES_BOTH"

MAIN_JOB_ID=$(sbatch --parsable \
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
        \"$SUMMARY_FILES_BOTH\" \
        \"$selected_simulations\" \
        \"$ANALYSIS_TYPE\" \
        \"$selected_rev_suffixes\" \
        \"$selected_models\" \
        \"$selected_configs_main\" \
        \"main\" \
        \"$COMBINED_PLOTS_DIR\"
    
    echo \"Main publication plots completed: \$(date)\"
    ")

# =============================================================================
# # Step 2: Supplementary Figure 1 (both simulations, 1:3 configs)
echo ""
echo "=== Step 2: Creating Supplementary Figure 1 (Both Simulations, 1:3 configs) ==="
echo "Target: Both simulations x 2 configs (1:3 sel and neu) = 4 panels"

selected_simulations_supp="tltt_08_20,gc_reentry_hunter"
selected_configs_supp1="config_ratio_1to3_sel,config_ratio_1to3_neu"
SUMMARY_FILE_SUPP="${BASE_PROJECT_ROOT}/tltt_08_20/results/${ANALYSIS_TYPE}/irrev/tree_analysis/all_results_summary.csv"

SUPP1_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${MAIN_JOB_ID} \
    --cpus-per-task=2 \
    --mem-per-cpu=6gb \
    --time=30:00 \
    --job-name=supp1-figure \
    --output="${LOG_FILE}" \
    --error="${LOG_FILE}" \
    --account=hoehnlab-share \
    --wrap="
    echo \"Supplementary figure 1 started: \$(date)\"
    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
    conda activate r_phylo
    cd \"$BASE_PROJECT_ROOT\"
    
    Rscript scripts/plotting/publication_plots_main.R \
        \"$SUMMARY_FILES_BOTH\" \
        \"$selected_simulations_supp\" \
        \"$ANALYSIS_TYPE\" \
        \"$selected_rev_suffixes\" \
        \"$selected_models\" \
        \"$selected_configs_supp1\" \
        \"supp_all_metrics\" \
        \"$COMBINED_PLOTS_DIR\"
    
    echo \"Supplementary figure 1 completed: \$(date)\"
    ")

# =============================================================================
# Step 3: Supplementary Figure 2 (both simulations, tree length analysis)
echo ""
echo "=== Step 3: Creating Supplementary Figure 2 (Tree Length Analysis) ==="
echo "Target: Both simulations x 4 configs (1:3 and 1:1, sel and neu) = 8 panels"

selected_configs_supp2="config_ratio_1to1_sel,config_ratio_1to1_neu,config_ratio_1to3_sel,config_ratio_1to3_neu"

echo "Summary files for supp2: $SUMMARY_FILES_BOTH"

SUPP2_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${SUPP1_JOB_ID} \
    --cpus-per-task=2 \
    --mem-per-cpu=6gb \
    --time=30:00 \
    --job-name=supp2-figure \
    --output="${LOG_FILE}" \
    --error="${LOG_FILE}" \
    --account=hoehnlab-share \
    --wrap="
    echo \"Supplementary figure 2 started: \$(date)\"
    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
    conda activate r_phylo
    cd \"$BASE_PROJECT_ROOT\"
    
    Rscript scripts/plotting/publication_plots_main.R \
        \"$SUMMARY_FILES_BOTH\" \
        \"$selected_simulations\" \
        \"$ANALYSIS_TYPE\" \
        \"$selected_rev_suffixes\" \
        \"$selected_models\" \
        \"$selected_configs_supp2\" \
        \"supp_tree_length\" \
        \"$COMBINED_PLOTS_DIR\"
    
    echo \"Supplementary figure 2 completed: \$(date)\"
    ")

echo ""
echo "=== Publication Plots Pipeline Submitted ==="
echo "Main figure job: $MAIN_JOB_ID"
echo "Supplementary figure 1 job: $SUPP1_JOB_ID" 
echo "Supplementary figure 2 job: $SUPP2_JOB_ID"
echo "Output directory: $COMBINED_PLOTS_DIR"
echo "Log file: $LOG_FILE"

# Final summary
echo ""
echo "=== Summary of Figures Being Created ==="
echo "1. Main Figure:"
echo "   - Simulations: tltt_08_20, gc_reentry_hunter"
echo "   - Configs: 1:1 sel, 1:1 neu (4 panels total)"
echo "   - Metrics: tree height, RF distance, MRCA accuracy"
echo "   - Facet labels: with simulation names"
echo ""
echo "2. Supplementary Figure 1:"
echo "   - Simulation: tltt_08_20, gc_reentry_hunter"
echo "   - Configs: 1:3 sel, 1:3 neu (4 panels total)"
echo "   - Metrics: tree height, RF distance, MRCA accuracy"
echo "   - Facet labels: with simulation names"
echo ""
echo "3. Supplementary Figure 2:"
echo "   - Simulations: tltt_08_20, gc_reentry_hunter"
echo "   - Configs: 1:1 sel, 1:1 neu, 1:3 sel, 1:3 neu (8 panels total)"
echo "   - Metric: tree length proportional error"
echo "   - Facet labels: with simulation names and config names"
echo ""
echo "=== Done ==="