#!/bin/bash
#SBATCH --job-name=simulate-data
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.out
#SBATCH --error=/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/slurm-output/%x/slurm-%A/slurm-%A_%a.err
#SBATCH --account=hoehnlab
#SBATCH --time=4:00:00
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

set -e

# Configuration
PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
simulation_run="tltt_08_20"
n_clones=20
n_processors=$n_clones
seed=12345

# Note: all analyses will use the same data as main_analysis

echo "=== Data Generation Configuration ==="
echo "Simulation run: $simulation_run"
echo "Number of clones: $n_clones"
echo "Seed: $seed"
echo "Note: sub_analysis and differentiation analysis will reuse main_analysis data"
echo "====================================="

# Directory setup
simulation_dir="${PROJECT_ROOT}/${simulation_run}"
log_dir="$simulation_dir/logs"
configs_dir="$simulation_dir/configs"
raw_data_dir="$simulation_dir/data/raw"

# Create necessary directories
mkdir -p "$simulation_dir" "$log_dir" "$configs_dir" "$raw_data_dir"

# Change to project directory
cd "$PROJECT_ROOT"

# Load conda environment
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate r_phylo

# Set up logging
log_file="$log_dir/01_simulation.log"
echo "=== Starting simulation at $(date) ===" | tee -a "$log_file"

# ------------------------------------------------------------------------
# Step 1: Generate configurations
echo "Step 1: Generating configurations..." | tee -a "$log_file"
# Generate configs
python scripts/beast/generate_simble_configs.py "$simulation_dir" >> "$log_file" 2>&1

# Count configs
config_count=$(find "$configs_dir" -name "*.json" 2>/dev/null | wc -l)
echo " Generated $config_count configuration files" | tee -a "$log_file"

if [[ $config_count -eq 0 ]]; then
    echo "ERROR: No configuration files generated!" | tee -a "$log_file"
    exit 1
fi

# ------------------------------------------------------------------------
# Step 2: Run simulations
echo "Step 2: Running simulations..." | tee -a "$log_file"

total_successful=0
total_failed=0

for config_file in "$configs_dir"/*.json; do
    config_name=$(basename "$config_file" .json)
    echo " Running: $config_name" | tee -a "$log_file"
    
    # Create data output directory
    output_dir="$raw_data_dir/$config_name"
    mkdir -p "$output_dir"
    
    # Build simble command
    simble_cmd="simble -o \"$output_dir\" -n $n_clones -p $n_processors --config \"$config_file\" --seed $seed --quiet"
    
    # Add migration rate if needed
    if [[ "$config_name" == *"GCmig_0.5"* ]]; then
        simble_cmd+=" --migration-rate 0.5"
    elif [[ "$config_name" == *"GCmig_"* ]]; then
        migration_rate=$(echo "$config_name" | grep -o 'GCmig_[0-9.]*' | cut -d'_' -f2)
        if [[ -n "$migration_rate" ]]; then
            simble_cmd+=" --migration-rate $migration_rate"
        fi
    fi
    
    # Run simulation
    if eval "$simble_cmd" >> "$log_file" 2>&1; then
        echo "    $config_name completed" | tee -a "$log_file"
        total_successful=$((total_successful + 1))
    else
        echo "    $config_name failed" | tee -a "$log_file"
        total_failed=$((total_failed + 1))
    fi
done

# Log summary
echo "" | tee -a "$log_file"
echo "=========== Summary ===========" | tee -a "$log_file"
echo "Configurations: $config_count" | tee -a "$log_file"
echo "Successful: $total_successful" | tee -a "$log_file"
echo "Failed: $total_failed" | tee -a "$log_file"
echo "===============================" | tee -a "$log_file"

if [[ $total_failed -gt 0 ]]; then
    echo "WARNING: $total_failed simulations failed" | tee -a "$log_file"
fi

echo "=== Simulations completed at $(date) ===" | tee -a "$log_file"