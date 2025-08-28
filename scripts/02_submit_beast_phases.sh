#!/bin/bash
# Configuration
PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
simulation_run="tltt_08_20"
analysis_scope="main_analysis"  # main_analysis, sub_analysis, differentiation_analysis
reversible="false"               # true or false
max_jobs_at_once=20              # how many jobs to run in parallel

# Script locations
run_script="${PROJECT_ROOT}/scripts/beast/run_beast_dowser.sh"
combo_script="${PROJECT_ROOT}/scripts/beast/create_beast_job_combinations.sh"

echo "=== Submitting BEAST Analysis Jobs ==="
echo "Simulation: $simulation_run"
echo "Analysis: $analysis_scope"
echo "Reversible: $reversible"
echo "======================================="

# Step 0: Generate job combinations
echo ""
echo "Step 0: Generating job combinations..."

# Determine raw data directory
raw_data_dir="${PROJECT_ROOT}/${simulation_run}/data/raw"
echo "Using raw data directory: $raw_data_dir"

if ! bash "$combo_script" "$raw_data_dir" "$simulation_run" "$analysis_scope"; then
    echo "✗ Failed to generate job combinations"
    exit 1
fi

echo "Job combinations generated successfully"

# Step 1: Submit GC Strict Clock Jobs (no dependencies)
# Count how many jobs we have
if [[ "$analysis_scope" == "main_analysis" ]]; then
    job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_gc_strict_clock.csv"
else
    job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_${analysis_scope}_gc_strict_clock.csv"
fi

if [[ ! -f "$job_file" ]]; then
    echo "ERROR: Job file not found: $job_file"
    exit 1
fi

gc_job_count=$(tail -n +2 "$job_file" | wc -l)
echo "Found $gc_job_count GC strict clock jobs"

# Submit the jobs
gc_job_id=$(sbatch --parsable --array=1-${gc_job_count}%${max_jobs_at_once} \
    $run_script "$simulation_run" "$analysis_scope" "gc_strict_clock" "$reversible")

if [[ $? -eq 0 ]]; then
    echo "GC strict clock jobs submitted with ID: $gc_job_id"
else
    echo "Failed to submit GC strict clock jobs"
    exit 1
fi

# Step 2: Submit TyCHE Models (depends on step 1)
echo ""
echo "Step 2: Submitting TyCHE model jobs (waiting for step 1)..."

# Count TyCHE jobs
if [[ "$analysis_scope" == "main_analysis" ]]; then
    tyche_job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_tyche_models.csv"
else
    tyche_job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_${analysis_scope}_tyche_models.csv"
fi

if [[ -f "$tyche_job_file" ]]; then
    tyche_job_count=$(tail -n +2 "$tyche_job_file" | wc -l)
    echo "Found $tyche_job_count TyCHE model jobs"
    
    # Submit with dependency on step 1
    tyche_job_id=$(sbatch --parsable --dependency=afterok:${gc_job_id} \
        --array=1-${tyche_job_count}%${max_jobs_at_once} \
        $run_script "$simulation_run" "$analysis_scope" "tyche_models" "$reversible")
    # sbatch --array=1-${tyche_job_count}%${max_jobs_at_once} \
    #     $run_script "$simulation_run" "$analysis_scope" "tyche_models" "$reversible"
    if [[ $? -eq 0 ]]; then
        echo "TyCHE model jobs submitted with ID: $tyche_job_id"
    else
        echo "Failed to submit TyCHE model jobs"
    fi
else
    echo "No TyCHE job file found, skipping..."
    tyche_job_id=""
fi

# Step 3: Submit Competing Models (depends on step 1)
# Count competing model jobs
if [[ "$analysis_scope" == "main_analysis" ]]; then
    competing_job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_competing_models.csv"
else
    competing_job_file="${PROJECT_ROOT}/${simulation_run}/configs/beast_job_combinations_${analysis_scope}_competing_models.csv"
fi

if [[ -f "$competing_job_file" ]]; then
    competing_job_count=$(tail -n +2 "$competing_job_file" | wc -l)
    echo "Found $competing_job_count competing model jobs"
    
    # Submit with dependency on step 1
    competing_job_id=$(sbatch --parsable --dependency=afterok:${gc_job_id} \
        --array=1-${competing_job_count}%${max_jobs_at_once} \
        $run_script "$simulation_run" "$analysis_scope" "competing_models" "$reversible")
    # sbatch --array=2-${competing_job_count} \
    #     $run_script "$simulation_run" "$analysis_scope" "competing_models" "$reversible"

    if [[ $? -eq 0 ]]; then
        echo "✓ Competing model jobs submitted with ID: $competing_job_id"
    else
        echo "✗ Failed to submit competing model jobs"
    fi
else
    echo "No competing models job file found, skipping..."
    competing_job_id=""
fi

# Summary
echo ""
echo "=== JOB SUBMISSION SUMMARY ==="
echo "GC Strict Clock: $gc_job_id ($gc_job_count jobs)"
[[ -n "$tyche_job_id" ]] && echo "TyCHE Models: $tyche_job_id ($tyche_job_count jobs)"
[[ -n "$competing_job_id" ]] && echo "Competing Models: $competing_job_id ($competing_job_count jobs)"
echo ""