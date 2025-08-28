#!/bin/bash
# scripts/beast/create_beast_job_combinations.sh

# Parse arguments
raw_data_dir="$1"
simulation_run="$2"
analysis_scope="$3"
output_dir="${4:-/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE/${simulation_run}/configs}"

if [[ ! -d "$raw_data_dir" ]]; then
    echo "ERROR: Raw data directory not found: $raw_data_dir"
    exit 1
fi

# Create output directory
mkdir -p "$output_dir"

# Analysis scope specific configurations
case "$analysis_scope" in
    main_analysis)
        echo "Configuring for MAIN ANALYSIS (all models, all phases)"
        
        # All templates for main analysis
        gc_strict_templates=("StrictClock_Standard")
        tyche_templates=(
            "ExpectedOccupancy_FixedClockRates"
            "ExpectedOccupancy_EstClockRates"
            "InstantSwitch_EstClockRates"
            "MixedSwitch_FixedClockRates"
            "MixedSwitch_EstClockRates"
        )
        competing_templates=(
            "StrictClock_AncestralReconstruction"
            "UCRelaxedClock_AncestralReconstruction"
        )
        germline_options=("true")
        file_suffix=""
        ;;
        
    sub_analysis)
        echo "Configuring for SUB ANALYSIS (time subsets)"
        
        # Reduced set for sub analysis
        gc_strict_templates=("StrictClock_Standard")
        tyche_templates=("ExpectedOccupancy_FixedClockRates")
        competing_templates=()
        germline_options=("true")
        file_suffix="_sub_analysis"
        # Define time subsets for sub analysis
        time_subsets=("50_100" "100_150" "150_200")
        ;;
        
    differentiation_analysis)
        echo "Configuring for DIFFERENTIATION ANALYSIS (3-state model only)"
        
        # Only the differentiation-specific model
        gc_strict_templates=()
        tyche_templates=("TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq_3state")     
        competing_templates=()
        germline_options=("true")
        file_suffix="_differentiation_analysis"
        ;;
        
    *)
        echo "ERROR: Invalid analysis_scope '$analysis_scope'"
        echo "Valid options: main_analysis, sub_analysis, differentiation_analysis"
        exit 1
        ;;
esac

# Define output files
gc_strict_file="${output_dir}/beast_job_combinations${file_suffix}_gc_strict_clock.csv"
tyche_file="${output_dir}/beast_job_combinations${file_suffix}_tyche_models.csv"
competing_file="${output_dir}/beast_job_combinations${file_suffix}_competing_models.csv"

# Create CSV headers with time_subset column for sub_analysis
if [[ "$analysis_scope" == "sub_analysis" ]]; then
    if [[ ${#gc_strict_templates[@]} -gt 0 ]]; then
        echo "config_name,template_id,germline,airr_file,time_subset" > "$gc_strict_file"
    fi
    if [[ ${#tyche_templates[@]} -gt 0 ]]; then
        echo "config_name,template_id,germline,airr_file,time_subset" > "$tyche_file"
    fi
else
    # Standard headers for other analyses
    if [[ ${#gc_strict_templates[@]} -gt 0 ]]; then
        echo "config_name,template_id,germline,airr_file" > "$gc_strict_file"
    fi
    if [[ ${#tyche_templates[@]} -gt 0 ]]; then
        echo "config_name,template_id,germline,airr_file" > "$tyche_file"
    fi
    if [[ ${#competing_templates[@]} -gt 0 ]]; then
        echo "config_name,template_id,germline,airr_file" > "$competing_file"
    fi
fi

echo "===== Generating BEAST Job Combinations ====="
echo "Raw data directory: $raw_data_dir"
echo "Simulation run: $simulation_run"
echo "Analysis scope: $analysis_scope"
echo "Output directory: $output_dir"
echo "=============================================="

# Initialize counters
gc_strict_count=0
tyche_count=0
competing_count=0

# Process each config directory
for config_dir in "$raw_data_dir"/config_*; do
    if [[ -d "$config_dir" ]]; then
        config_name=$(basename "$config_dir")

        # Filter configs based on analysis scope
        case "$analysis_scope" in
            main_analysis)
                # Process all configs
                ;;
            sub_analysis)
                # Only process sel and neu configs
                if [[ ! "$config_name" =~ _1to3_(sel|neu)$ ]]; then
                    echo "  Skipping $config_name (not sel/neu for sub_analysis)"
                    continue
                fi
                ;;
            differentiation_analysis)
                # Only process sel configs
                if [[ ! "$config_name" =~ _1to3_sel$ ]]; then
                    echo "  Skipping $config_name (not sel for differentiation_analysis)"
                    continue
                fi
                ;;
        esac

        airr_file="$config_dir/all_samples_airr.tsv"
        
        echo "Processing config: $config_name"
        
        # Check if the AIRR file exists
        if [[ -f "$airr_file" ]]; then
            echo "  ✓ Found AIRR file: $airr_file"
            
            if [[ "$analysis_scope" == "sub_analysis" ]]; then
                echo "  Creating jobs for time subsets: ${time_subsets[*]}"
                
                # For sub_analysis, create jobs for each time subset
                for time_subset in "${time_subsets[@]}"; do
                    # GC STRICT CLOCK JOBS
                    if [[ ${#gc_strict_templates[@]} -gt 0 ]]; then
                        for template in "${gc_strict_templates[@]}"; do
                            for germline in "${germline_options[@]}"; do
                                echo "$config_name,$template,$germline,$airr_file,$time_subset" >> "$gc_strict_file"
                                ((gc_strict_count++))
                            done
                        done
                    fi
                    
                    # TYCHE MODEL JOBS
                    if [[ ${#tyche_templates[@]} -gt 0 ]]; then
                        for template in "${tyche_templates[@]}"; do
                            for germline in "${germline_options[@]}"; do
                                echo "$config_name,$template,$germline,$airr_file,$time_subset" >> "$tyche_file"
                                ((tyche_count++))
                            done
                        done
                    fi
                done
            else
                # For main_analysis and differentiation_analysis, standard job creation
                
                # GC STRICT CLOCK JOBS
                if [[ ${#gc_strict_templates[@]} -gt 0 ]]; then
                    for template in "${gc_strict_templates[@]}"; do
                        for germline in "${germline_options[@]}"; do
                            echo "$config_name,$template,$germline,$airr_file" >> "$gc_strict_file"
                            ((gc_strict_count++))
                        done
                    done
                fi
                
                # TYCHE MODEL JOBS
                if [[ ${#tyche_templates[@]} -gt 0 ]]; then
                    for template in "${tyche_templates[@]}"; do
                        for germline in "${germline_options[@]}"; do
                            echo "$config_name,$template,$germline,$airr_file" >> "$tyche_file"
                            ((tyche_count++))
                        done
                    done
                fi
                
                # COMPETING MODEL JOBS
                if [[ ${#competing_templates[@]} -gt 0 ]]; then
                    for template in "${competing_templates[@]}"; do
                        for germline in "${germline_options[@]}"; do
                            echo "$config_name,$template,$germline,$airr_file" >> "$competing_file"
                            ((competing_count++))
                        done
                    done
                fi
            fi
            
        else
            echo "  WARNING: AIRR file not found: $airr_file"
        fi
    fi
done

# Summary
# Count actual jobs in files (excluding headers)
if [[ -f "$gc_strict_file" ]]; then
    gc_strict_jobs=$(tail -n +2 "$gc_strict_file" | wc -l)
else
    gc_strict_jobs=0
fi

if [[ -f "$tyche_file" ]]; then
    tyche_jobs=$(tail -n +2 "$tyche_file" | wc -l)
else
    tyche_jobs=0
fi

if [[ -f "$competing_file" ]]; then
    competing_jobs=$(tail -n +2 "$competing_file" | wc -l)
else
    competing_jobs=0
fi

total_jobs=$((gc_strict_jobs + tyche_jobs + competing_jobs))

echo ""
echo "=== Job Generation Summary ==="
if [[ ${#gc_strict_templates[@]} -gt 0 ]]; then
    echo "GC Strict Clock: $gc_strict_count jobs → $(basename "$gc_strict_file")"
fi
if [[ ${#tyche_templates[@]} -gt 0 ]]; then
    echo "TyCHE Models: $tyche_count jobs → $(basename "$tyche_file")"
fi
if [[ ${#competing_templates[@]} -gt 0 ]]; then
    echo "Competing Models: $competing_count jobs → $(basename "$competing_file")"
fi

echo "Total jobs: $total_jobs"

# Validate job counts
if [[ $total_jobs -eq 0 ]]; then
    echo ""
    echo "ERROR: No jobs generated!"
    exit 1
fi

echo ""
echo "=== Template Summary ==="
echo "GC Strict Clock templates: ${gc_strict_templates[*]:-none}"
echo "TyCHE model templates: ${tyche_templates[*]:-none}"
echo "Competing model templates: ${competing_templates[*]:-none}"

echo ""
echo "✓ Job combinations generated successfully!"