#!/bin/bash
# master_pipeline.sh

SIMULATION_NAME="${1:-tltt_08_20}"

echo "=== Starting Master Pipeline: $(date) ==="
echo "Simulation: $SIMULATION_NAME"

# Step 1: Simulate data
echo ""
echo "=== Step 1: Data Simulation ==="
sbatch scripts/01_simulate_data.sh "$SIMULATION_NAME"

# Step 2: Submit BEAST jobs in multiple phases for GC strict clock and TyCHE/competing models
echo ""
echo "=== Step 2: BEAST Analysis ==="
bash scripts/02_submit_beast_phases.sh "$SIMULATION_NAME" main_analysis irrev
bash scripts/02_submit_beast_phases.sh "$SIMULATION_NAME" main_analysis rev
bash scripts/02_submit_beast_phases.sh "$SIMULATION_NAME" differentiation_analysis irrev

# Step 3-1: Tree analysis for main_analysis
echo ""
echo "=== Step 3-1: Tree Analysis (main_analysis) ==="
bash scripts/03_analyze_trees.sh "$SIMULATION_NAME" irrev main_analysis
bash scripts/03_analyze_trees.sh "$SIMULATION_NAME" rev main_analysis

# Step 3-2: Tree analysis for GC re-entry simulation
echo ""
echo "=== Step 3-2: Tree Analysis (gc_reentry_hunter) ==="
bash scripts/setup/import_gc_reentry_data.sh config_ratio_1to1_sel 8_28
bash scripts/setup/import_gc_reentry_data.sh config_ratio_1to1_neu uniform_neutral_8_29
bash scripts/03_analyze_trees.sh gc_reentry_hunter rev main_analysis

# Step 3-3: Tree analysis for differentiation_analysis
echo ""
echo "=== Step 3-3: Tree Analysis (differentiation_analysis) ==="
bash scripts/03_analyze_trees.sh "$SIMULATION_NAME" irrev differentiation_analysis

# Step 4: Wait for all analyses to complete, then create publication plots
bash scripts/04_generate_figures.sh main_analysis
bash scripts/04_generate_figures.sh differentiation_analysis

echo ""
echo "=== Master Pipeline Completed: $(date) ==="