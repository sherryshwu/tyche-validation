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

# Step 3: Tree analysis for main_analysis
echo ""
echo "=== Step 3: Tree Analysis (main_analysis) ==="
bash scripts/03_tree_analysis_main.sh "$SIMULATION_NAME" irrev
bash scripts/03_tree_analysis_main.sh "$SIMULATION_NAME" rev

# Step 4: Tree analysis for GC re-entry simulation
echo ""
echo "=== Step 4: Tree Analysis (gc_reentry_hunter) ==="
bash scripts/03_tree_analysis_main.sh gc_reentry_hunter irrev

# Step 5: Tree analysis for differentiation_analysis
echo ""
echo "=== Step 5: Tree Analysis (differentiation_analysis) ==="
bash scripts/03_tree_analysis_differentiation.sh "$SIMULATION_NAME" irrev

# To be finalized...
echo "=== Step 6: Tree Analysis (sub_analysis) ==="
bash scripts/03_tree_analysis_sub.sh "$SIMULATION_NAME" irrev

echo ""
echo "=== Master Pipeline Completed: $(date) ==="