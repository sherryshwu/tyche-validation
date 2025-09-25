#!/bin/bash
# setup_gc_reentry.sh
CONFIG_NAME="${1:-config_ratio_1to1_sel}"
DIR_SUFFIX="${2:-8_28}"

PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
HUNTER_SIMULATION_NAME="gc_reentry_hunter"
REV_SUFFIX="rev"
ANALYSIS_TYPE="main_analysis"
config_name="$CONFIG_NAME"

# Hunter's source directories
HUNTER_SIMBLE_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/simble_sims_gc_reentry_${DIR_SUFFIX}/"
HUNTER_BEAST_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/beast_results_${DIR_SUFFIX}/"

# Conditional dowser source path based on config name
if [[ "$config_name" == *"1to1_sel"* ]]; then
    HUNTER_DOWSER_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/output_${DIR_SUFFIX}_sherry/"
    echo "Using Sherry-specific dowser source for config: $config_name"
else
    HUNTER_DOWSER_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/output_${DIR_SUFFIX}/"
    echo "Using standard dowser source for config: $config_name"
fi

# Setup Hunter's directory structure
HUNTER_BASE_DATA_DIR="${PROJECT_ROOT}/${HUNTER_SIMULATION_NAME}/data"
HUNTER_RAW_DATA_DIR="${HUNTER_BASE_DATA_DIR}/raw/${config_name}"
HUNTER_PROCESSED_DATA_DIR="${HUNTER_BASE_DATA_DIR}/processed"
HUNTER_RESULTS_BASE_DIR="${PROJECT_ROOT}/${HUNTER_SIMULATION_NAME}/results/${ANALYSIS_TYPE}/${REV_SUFFIX}"
HUNTER_DOWSER_DIR="${HUNTER_PROCESSED_DATA_DIR}/${ANALYSIS_TYPE}/dowser_processed_trees/${REV_SUFFIX}"

echo "=== Setting up Hunter's GC Reentry Analysis ==="
echo "Hunter simulation: $HUNTER_SIMULATION_NAME"

# 1. Setup true trees
echo "Setting up true trees..."
mkdir -p "$HUNTER_RAW_DATA_DIR"
if [[ -d "$HUNTER_SIMBLE_SOURCE" ]]; then
    rsync -av "$HUNTER_SIMBLE_SOURCE" "$HUNTER_RAW_DATA_DIR/"
    echo "✓ True trees copied to: $HUNTER_RAW_DATA_DIR"
fi

# 2. Setup dowser trees
echo "Setting up dowser trees..."
mkdir -p "$HUNTER_DOWSER_DIR"
if [[ -d "$HUNTER_DOWSER_SOURCE" ]]; then
    echo "Copying and renaming dowser files..."
    
    # Copy all files first
    rsync -av "$HUNTER_DOWSER_SOURCE" "${HUNTER_DOWSER_DIR}/temp_dowser/"
    
    # Process each file in the temp directory
    find "${HUNTER_DOWSER_DIR}/temp_dowser/" -type f -name "*.rds" | while read file; do
        filename=$(basename "$file")
        
        # 1. Capitalize first letters
        new_filename=$(echo "$filename" | sed 's/expectedOccupancy/ExpectedOccupancy/' | sed 's/expectedoccupancy/ExpectedOccupancy/' | sed 's/instantSwitch/InstantSwitch/' | sed 's/instantswitch/InstantSwitch/' | sed 's/mixedSwitch/MixedSwitch/' | sed 's/mixedswitch/MixedSwitch/' | sed 's/strictClock/StrictClock/' | sed 's/strictclock/StrictClock/' | sed 's/ucRelaxedClock/UCRelaxedClock/' | sed 's/ucrelaxedclock/UCRelaxedClock/')
        
        # 2. Replace EstTraitClockRates with EstClockRates and FixedTraitClockRates with FixedClockRates
        new_filename=$(echo "$new_filename" | sed 's/EstTraitClockRates/EstClockRates/g' | sed 's/FixedTraitClockRates/FixedClockRates/g')
        
        # 3. Remove _EmpFreq
        new_filename=$(echo "$new_filename" | sed 's/_EmpFreq//g')
        
        # 4. Replace _1_2_3_4_5_6..._20 with _all (matches any sequence of numbers with underscores)
        new_filename=$(echo "$new_filename" | sed 's/_[0-9_]*[0-9]_trees\.rds$/_all_trees.rds/')

        # 5. Remove gc_reentry_sims_
        new_filename=$(echo "$new_filename" | sed 's/gc_reentry_sims_//g')

        # 6. Remove uniform_neutral
        new_filename=$(echo "$new_filename" | sed 's/uniform_neutral//g')

        # 6. Add config_name to the beginning
        new_filename="${config_name}_${new_filename}"
        
        # Copy with new name
        cp "$file" "${HUNTER_DOWSER_DIR}/${new_filename}"
    done
    
    # Remove temp directory
    rm -rf "${HUNTER_DOWSER_DIR}/temp_dowser/"
    
    echo "✓ Dowser trees copied and renamed to: $HUNTER_DOWSER_DIR"
else
    echo "⚠ Dowser source not found (may still be running): $HUNTER_DOWSER_SOURCE"
fi

# 3. Setup BEAST results with model type inference
echo "Setting up BEAST results..."
mkdir -p "${HUNTER_RESULTS_BASE_DIR}/tyche_models/beast_raw_output/${config_name}"
mkdir -p "${HUNTER_RESULTS_BASE_DIR}/competing_models/beast_raw_output/${config_name}"

if [[ -d "$HUNTER_BEAST_SOURCE" ]]; then
    echo "Copying all BEAST files and organizing by model type..."
    
    # Copy tyche_models files (case-insensitive)
    find "$HUNTER_BEAST_SOURCE" -iname "ExpectedOccupancy*" -o -iname "InstantSwitch*" -o -iname "MixedSwitch*" | while read file; do
        cp "$file" "${HUNTER_RESULTS_BASE_DIR}/tyche_models/beast_raw_output/${config_name}/"
    done
    
    # Copy competing_models files (case-insensitive)  
    find "$HUNTER_BEAST_SOURCE" -iname "strictclock*" -o -iname "ucrelaxedclock*" | while read file; do
        cp "$file" "${HUNTER_RESULTS_BASE_DIR}/competing_models/beast_raw_output/${config_name}/"
    done
    
    echo "BEAST results organized by model type"

    # Process tyche_models files
    cd "${HUNTER_RESULTS_BASE_DIR}/tyche_models/beast_raw_output/${config_name}"
    for file in *; do
        if [[ -f "$file" ]]; then
            new_name="$file"
            
            # 1. Capitalize first letters
            new_name=$(echo "$new_name" | sed 's/^expectedOccupancy/ExpectedOccupancy/' | sed 's/^expectedoccupancy/ExpectedOccupancy/' | sed 's/^instantSwitch/InstantSwitch/' | sed 's/^instantswitch/InstantSwitch/' | sed 's/^mixedSwitch/MixedSwitch/' | sed 's/^mixedswitch/MixedSwitch/')
            
            # 2. Replace EstTraitClockRates with EstClockRates and FixedTraitClockRates with FixedClockRates
            new_name=$(echo "$new_name" | sed 's/EstTraitClockRates/EstClockRates/g' | sed 's/FixedTraitClockRates/FixedClockRates/g')
            
            # 3. Remove _EmpFreq
            new_name=$(echo "$new_name" | sed 's/_EmpFreq//g')
            
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        fi
    done
    
    # Process competing_models files
    cd "${HUNTER_RESULTS_BASE_DIR}/competing_models/beast_raw_output/${config_name}"
    for file in *; do
        if [[ -f "$file" ]]; then
            new_name="$file"
            
            # 1. Capitalize first letters
            new_name=$(echo "$new_name" | sed 's/^strictClock/StrictClock/' | sed 's/^strictclock/StrictClock/' | sed 's/^ucRelaxedClock/UCRelaxedClock/' | sed 's/^ucrelaxedclock/UCRelaxedClock/')
            
            # 2. Replace EstTraitClockRates with EstClockRates and FixedTraitClockRates with FixedClockRates
            new_name=$(echo "$new_name" | sed 's/EstTraitClockRates/EstClockRates/g' | sed 's/FixedTraitClockRates/FixedClockRates/g')
            
            # 3. Remove _EmpFreq
            new_name=$(echo "$new_name" | sed 's/_EmpFreq//g')
            
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        fi
    done

    echo "BEAST results organized by model type and filenames standardized"
fi

# 4. Create configs directory
HUNTER_CONFIGS_DIR="${PROJECT_ROOT}/${HUNTER_SIMULATION_NAME}/configs"
mkdir -p "$HUNTER_CONFIGS_DIR"

echo ""
echo "✓ Hunter's data setup complete!"
echo ""
echo "Directory structure created:"
echo "  Data: $HUNTER_BASE_DATA_DIR"
echo "  Results: $HUNTER_RESULTS_BASE_DIR"
echo "  Configs: $HUNTER_CONFIGS_DIR"
echo ""
echo "File naming changes applied:"
echo "  - EstTraitClockRates → EstClockRates"
echo "  - FixedTraitClockRates → FixedClockRates"
echo "  - Removed _EmpFreq suffix"
echo "  - Capitalized first letters"
echo "  - Dowser files: added config prefix and _1_2_3...20 → _all"