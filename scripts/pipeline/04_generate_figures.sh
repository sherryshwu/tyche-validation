#!/bin/bash
# Multi-analysis publication plotting script

# Get parameters
ANALYSIS_TYPE="${1:-main_analysis}"
SIM_INPUT="${2:-tltt_12_19,gc_reentry_12_18}"

BASE_PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"

SELECTED_SIMULATIONS_LIST=$(echo $SIM_INPUT | tr ',' ' ')

echo "=== Creating Publication Plots: $(date) ==="
echo "Analysis type: $ANALYSIS_TYPE"

# Setup combined output directory
COMBINED_PLOTS_DIR="${BASE_PROJECT_ROOT}/figures/dec2025/${ANALYSIS_TYPE}"
LOGS_DIR="${BASE_PROJECT_ROOT}/logs/dec2025/publication_plots"
mkdir -p "$COMBINED_PLOTS_DIR" "$LOGS_DIR"

# Setup logging
MAIN_LOG_FILE="${LOGS_DIR}/create_main_figures_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$MAIN_LOG_FILE")
exec 2>&1

# Common settings
selected_models="EO_Est,SC_AR,UCLD_AR"
selected_rev_suffixes="irrev,rev"

case "$ANALYSIS_TYPE" in
    "main_analysis")
        echo "=== Creating Main Analysis Publication Plots ==="
        # =============================================================================
        # Step 1: Main Publication Plots (Cross-simulation comparison)
        echo ""
        echo "=== Step 1: Creating Main Publication Plots ==="
        echo "Target: Two simulations x One 1:1 config = Two panels"

        selected_configs_main="config_ratio_1to1_sel"
        selected_configs_supp="config_ratio_1to1_neu"

        # Build summary file paths for both simulations
        SUMMARY_FILES_BOTH=""

        for sim in $SELECTED_SIMULATIONS_LIST; do
            if [[ "$sim" == *"gc_reentry"* ]]; then rev_suffix="rev"; else rev_suffix="irrev"; fi
            summary_file="${BASE_PROJECT_ROOT}/${sim}/results/${ANALYSIS_TYPE}/${rev_suffix}/tree_analysis/all_results_summary.csv"
            if [[ -f "$summary_file" ]]; then
                if [[ -n "$SUMMARY_FILES_BOTH" ]]; then
                    SUMMARY_FILES_BOTH="${SUMMARY_FILES_BOTH},${summary_file}"
                else
                    SUMMARY_FILES_BOTH="$summary_file"
                fi
            fi
        done

        echo "Summary files for main: $SUMMARY_FILES_BOTH"

        # Create separate log files for each step
        STEP1_LOG="${LOGS_DIR}/step1_main_plots_$(date +%Y%m%d_%H%M%S).log"
        STEP2_LOG="${LOGS_DIR}/step2_supp_plots_$(date +%Y%m%d_%H%M%S).log"

        MAIN_JOB_ID=$(sbatch --parsable \
            --cpus-per-task=2 \
            --mem-per-cpu=6gb \
            --time=45:00 \
            --job-name=main-fig-1to1-sel \
            --output="${STEP1_LOG}" \
            --error="${STEP1_LOG}" \
            --account=hoehnlab-share \
            --wrap="
            echo \"Main publication plots started: \$(date)\"
            source /optnfs/common/miniconda3/etc/profile.d/conda.sh
            conda activate r_phylo_4.4
            cd \"$BASE_PROJECT_ROOT\"
            
            Rscript scripts/visualization/create_main_figures.R \
                \"$SUMMARY_FILES_BOTH\" \
                \"$SIM_INPUT\" \
                \"$ANALYSIS_TYPE\" \
                \"$selected_rev_suffixes\" \
                \"$selected_models\" \
                \"$selected_configs_main\" \
                \"main\" \
                \"$COMBINED_PLOTS_DIR\"

            echo \"Main publication plots completed: \$(date)\"
            ")

        # =============================================================================
        # Step 2: Supplementary Figure 1 (both simulations, 1:1 uniform neutral)
        echo ""
        echo "=== Step 2: Creating Supplementary Figure 1 (Uniform neutral in both simulations) ==="
        echo "Target: Both simulations x 1 configs (1:1 neu) = 2 panels"

        SUMMARY_FILE_SUPP="${BASE_PROJECT_ROOT}/tltt_12_19/results/${ANALYSIS_TYPE}/irrev/tree_analysis/all_results_summary.csv"

        SUPP1_JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${MAIN_JOB_ID} \
            --cpus-per-task=2 \
            --mem-per-cpu=6gb \
            --time=30:00 \
            --job-name=supp-fig-1to1-neu \
            --output="${STEP2_LOG}" \
            --error="${STEP2_LOG}" \
            --account=hoehnlab-share \
            --wrap="
            echo \"Supplementary figure started: \$(date)\"
            source /optnfs/common/miniconda3/etc/profile.d/conda.sh
            conda activate r_phylo_4.4
            cd \"$BASE_PROJECT_ROOT\"
            
            Rscript scripts/visualization/create_main_figures.R \
                \"$SUMMARY_FILES_BOTH\" \
                \"$SIM_INPUT\" \
                \"$ANALYSIS_TYPE\" \
                \"$selected_rev_suffixes\" \
                \"$selected_models\" \
                \"$selected_configs_supp\" \
                \"supp\" \
                \"$COMBINED_PLOTS_DIR\"
            
            echo \"Supplementary figure completed: \$(date)\"
            ")

        echo ""
        echo "=== Main Analysis Publication Plots Submitted ==="
        echo "Main figure (1:1 Selection) job: $MAIN_JOB_ID"
        echo "Supplementary figure (1:1 Neutral) job: $SUPP1_JOB_ID" 
        echo "Output directory: $COMBINED_PLOTS_DIR"

        echo ""
        echo "=== Summary of Log Files ==="
        echo "Main script log: $MAIN_LOG_FILE"
        echo "Step 1 (main plots): $STEP1_LOG"
        echo "Step 2 (supp figure 1): $STEP2_LOG"

        echo ""
        echo "=== Summary of Figures Being Created ==="
        echo "1. Main Figure:"
        echo "   - Simulations: tltt_08_20, gc_reentry"
        echo "   - Configs: 1:1 sel (8 panels total)"
        echo "   - Metrics: tree height, tree length, RF distance, MRCA accuracy"
        echo "   - Facet labels: with simulation names"
        echo ""
        echo "2. Supplementary Figure 1:"
        echo "   - Simulation: tltt_08_20"
        echo "   - Configs: 1:1 neu (8 panels total)"
        echo "   - Metrics: tree height, tree length, RF distance, MRCA accuracy"
        echo "   - Facet labels: with simulation names"
        echo ""
        ;;

    "differentiation_analysis") 
        echo "=== Creating Differentiation Analysis Publication Plots ==="
        # =============================================================================
        # Step 4: Create Differentiation Timing Plots (for differentiation_analysis)
        echo ""
        echo "=== Step 3: Creating Differentiation Timing Plots ==="
        
        # Find the job list file for differentiation analysis
        JOB_LIST_FILE="${BASE_PROJECT_ROOT}/tltt_08_20/configs/tree_analysis_jobs_differentiation_analysis_irrev.csv"
        
        if [[ -f "$JOB_LIST_FILE" ]]; then
            # Read the total jobs from the file
            TOTAL_JOBS=$(tail -n +2 "$JOB_LIST_FILE" | wc -l)
            
            DIFF_TIMING_LOG_DIR="${BASE_PROJECT_ROOT}/tltt_08_20/logs/differentiation_analysis/differentiation_timing_jobs"
            mkdir -p "$DIFF_TIMING_LOG_DIR"
            
            DIFF_TIMING_JOB_ID=$(sbatch --parsable \
                --array=1-${TOTAL_JOBS}%5 \
                --cpus-per-task=2 \
                --mem-per-cpu=4gb \
                --time=30:00 \
                --job-name=diff-timing-differentiation_analysis \
                --output="${DIFF_TIMING_LOG_DIR}/diff_timing_%A_%a.out" \
                --error="${DIFF_TIMING_LOG_DIR}/diff_timing_%A_%a.err" \
                --account=hoehnlab-share \
                --wrap="
                    echo \"Differentiation timing job \$SLURM_ARRAY_TASK_ID started: \$(date)\"
                    source /optnfs/common/miniconda3/etc/profile.d/conda.sh
                    conda activate r_phylo_4.4
                    cd \"$BASE_PROJECT_ROOT\"
                    
                    Rscript scripts/visualization/create_differentiation_figures.R \
                        \"$JOB_LIST_FILE\" \
                        \"\$SLURM_ARRAY_TASK_ID\" \
                        \"${BASE_PROJECT_ROOT}/tltt_08_20/results/differentiation_analysis/irrev/tree_analysis\" \
                        \"$COMBINED_PLOTS_DIR\"
                    
                    echo \"Differentiation timing job \$SLURM_ARRAY_TASK_ID completed: \$(date)\"
                ")
            
            echo "Differentiation timing jobs submitted with ID: $DIFF_TIMING_JOB_ID"
            echo "Job logs in: $DIFF_TIMING_LOG_DIR"
        else
            echo "Warning: Job list file not found for differentiation timing plots: $JOB_LIST_FILE"
        fi

        echo ""
        echo "=== Differentiation Analysis Publication Plots Submitted ==="
        echo "Differentiation timing job: $DIFF_TIMING_JOB_ID"
        echo "Output directory: $COMBINED_PLOTS_DIR"

        echo ""
        echo "=== Summary of Figures Being Created ==="
        echo "1. Main Figure for Differentiation Timing Analysis:"
        echo "   - Simulations: tltt_08_20"
        echo "   - Configs: 1:3 sel (2 panels total)"
        ;;
    
    *)
        echo "Error: Unknown analysis type '$ANALYSIS_TYPE'"
        echo "Usage: $0 [main_analysis|differentiation_analysis]"
        exit 1
        ;;
esac

echo ""
echo "=== Publication Plots Pipeline Submitted ==="
echo "Log file: $MAIN_LOG_FILE"
echo "=== Done ==="