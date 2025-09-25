#!/bin/bash
# Multi-analysis publication plotting script

# Get parameters
ANALYSIS_TYPE="${1:-main_analysis}"
BASE_PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"

echo "=== Creating Publication Plots: $(date) ==="
echo "Analysis type: $ANALYSIS_TYPE"

# Setup combined output directory
COMBINED_PLOTS_DIR="${BASE_PROJECT_ROOT}/figures/${ANALYSIS_TYPE}"
LOGS_DIR="${BASE_PROJECT_ROOT}/logs/publication_plots"
mkdir -p "$COMBINED_PLOTS_DIR" "$LOGS_DIR"

# Setup logging
MAIN_LOG_FILE="${LOGS_DIR}/create_main_figures_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$MAIN_LOG_FILE")
exec 2>&1

# Common settings
selected_models="EO_Fixed,EO_Est,IS_Est,SC_AR,UCLD_AR"
selected_rev_suffixes="irrev,rev"

case "$ANALYSIS_TYPE" in
    "main_analysis")
        echo "=== Creating Main Analysis Publication Plots ==="
        # =============================================================================
        # Step 1: Main Publication Plots (Cross-simulation comparison)
        echo ""
        echo "=== Step 1: Creating Main Publication Plots ==="
        echo "Target: Two simulations x Two 1:1 configs = Four panels"

        selected_configs_main="config_ratio_1to1_sel,config_ratio_1to1_neu"
        selected_simulations="tltt_08_20,gc_reentry_hunter"

        # Build summary file paths for both simulations
        SUMMARY_FILES_BOTH=""
        declare -A sim_rev_map
        sim_rev_map[tltt_08_20]="irrev"
        sim_rev_map[gc_reentry_hunter]="rev"

        for sim in tltt_08_20 gc_reentry_hunter; do
            rev_suffix=${sim_rev_map[$sim]}
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
        STEP2_LOG="${LOGS_DIR}/step2_supp1_plots_$(date +%Y%m%d_%H%M%S).log"
        STEP3_LOG="${LOGS_DIR}/step3_supp2_plots_$(date +%Y%m%d_%H%M%S).log"

        MAIN_JOB_ID=$(sbatch --parsable \
            --cpus-per-task=2 \
            --mem-per-cpu=6gb \
            --time=45:00 \
            --job-name=main-publication-plots \
            --output="${STEP1_LOG}" \
            --error="${STEP1_LOG}" \
            --account=hoehnlab-share \
            --wrap="
            echo \"Main publication plots started: \$(date)\"
            source /optnfs/common/miniconda3/etc/profile.d/conda.sh
            conda activate r_phylo
            cd \"$BASE_PROJECT_ROOT\"
            
            Rscript scripts/visualization/create_main_figures.R \
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
        # Step 2: Supplementary Figure 1 (both simulations, 1:3 configs)
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
            --output="${STEP2_LOG}" \
            --error="${STEP2_LOG}" \
            --account=hoehnlab-share \
            --wrap="
            echo \"Supplementary figure 1 started: \$(date)\"
            source /optnfs/common/miniconda3/etc/profile.d/conda.sh
            conda activate r_phylo
            cd \"$BASE_PROJECT_ROOT\"
            
            Rscript scripts/visualization/create_main_figures.R \
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
            --output="${STEP3_LOG}" \
            --error="${STEP3_LOG}" \
            --account=hoehnlab-share \
            --wrap="
            echo \"Supplementary figure 2 started: \$(date)\"
            source /optnfs/common/miniconda3/etc/profile.d/conda.sh
            conda activate r_phylo
            cd \"$BASE_PROJECT_ROOT\"
            
            Rscript scripts/visualization/create_main_figures.R \
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
        echo "=== Main Analysis Publication Plots Submitted ==="
        echo "Main figure job: $MAIN_JOB_ID"
        echo "Supplementary figure 1 job: $SUPP1_JOB_ID" 
        echo "Supplementary figure 2 job: $SUPP2_JOB_ID"
        echo "Output directory: $COMBINED_PLOTS_DIR"

        echo ""
        echo "=== Summary of Log Files ==="
        echo "Main script log: $MAIN_LOG_FILE"
        echo "Step 1 (main plots): $STEP1_LOG"
        echo "Step 2 (supp figure 1): $STEP2_LOG"
        echo "Step 3 (supp figure 2): $STEP3_LOG"

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
        ;;

    "differentiation_analysis") 
        echo "=== Creating Differentiation Analysis Publication Plots ==="
        # =============================================================================
        # Step 4: Create Differentiation Timing Plots (for differentiation_analysis)
        echo ""
        echo "=== Step 4: Creating Differentiation Timing Plots ==="
        
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
                    conda activate r_phylo
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