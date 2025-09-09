#!/bin/bash
# setup_gc_reentry.sh
CONFIG_NAME="${1:-config_ratio_1to1_sel}"
DIR_SUFFIX="${2:-8_28}"

PROJECT_ROOT="/dartfs/rc/lab/H/HoehnK/Sherry/beast_workspace/TyCHE"
HUNTER_SIMULATION_NAME="gc_reentry_hunter"
REV_SUFFIX="irrev"
ANALYSIS_TYPE="main_analysis"
config_name="$CONFIG_NAME"

# Hunter's source directories
HUNTER_SIMBLE_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/simble_sims_gc_reentry_${DIR_SUFFIX}/"
HUNTER_BEAST_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/beast_results_${DIR_SUFFIX}/"
HUNTER_DOWSER_SOURCE="/dartfs/rc/lab/H/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/output_${DIR_SUFFIX}_sherry/"

# Setup Hunter's directory structure
HUNTER_BASE_DATA_DIR="${PROJECT_ROOT}/${HUNTER_SIMULATION_NAME}/data"
HUNTER_RAW_DATA_DIR="${HUNTER_BASE_DATA_DIR}/raw/${config_name}"
HUNTER_PROCESSED_DATA_DIR="${HUNTER_BASE_DATA_DIR}/processed"
HUNTER_RESULTS_BASE_DIR="${PROJECT_ROOT}/${HUNTER_SIMULATION_NAME}/results/${ANALYSIS_TYPE}/${REV_SUFFIX}"

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
HUNTER_DOWSER_DIR="${HUNTER_PROCESSED_DATA_DIR}/${ANALYSIS_TYPE}/dowser_processed_trees/${REV_SUFFIX}"
mkdir -p "$HUNTER_DOWSER_DIR"
if [[ -d "$HUNTER_DOWSER_SOURCE" ]]; then
    rsync -av "$HUNTER_DOWSER_SOURCE" "$HUNTER_DOWSER_DIR/"
    echo "✓ Dowser trees copied to: $HUNTER_DOWSER_DIR"
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

    # Capitalize first letters of filenames in tyche_models
    cd "${HUNTER_RESULTS_BASE_DIR}/tyche_models/beast_raw_output/${config_name}"
    for file in expectedOccupancy* expectedoccupancy* instantSwitch* instantswitch* mixedSwitch* mixedswitch*; do
        if [[ -f "$file" ]]; then
            new_name=$(echo "$file" | sed 's/^expectedOccupancy/ExpectedOccupancy/' | sed 's/^expectedoccupancy/ExpectedOccupancy/' | sed 's/^instantSwitch/InstantSwitch/' | sed 's/^instantswitch/InstantSwitch/' | sed 's/^mixedSwitch/MixedSwitch/' | sed 's/^mixedswitch/MixedSwitch/')
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        fi
    done
    
    # Capitalize first letters of filenames in competing_models
    cd "${HUNTER_RESULTS_BASE_DIR}/competing_models/beast_raw_output/${config_name}"
    for file in strictClock* strictclock* ucRelaxedClock* ucrelaxedclock*; do
        if [[ -f "$file" ]]; then
            new_name=$(echo "$file" | sed 's/^strictClock/StrictClock/' | sed 's/^strictclock/StrictClock/' | sed 's/^ucRelaxedClock/UCRelaxedClock/' | sed 's/^ucrelaxedclock/UCRelaxedClock/')
            if [[ "$file" != "$new_name" ]]; then
                mv "$file" "$new_name"
            fi
        fi
    done

    echo "BEAST results organized by model type and filenames capitalized"
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