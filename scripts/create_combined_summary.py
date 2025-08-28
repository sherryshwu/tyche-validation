#!/usr/bin/env python3

# Create combined summary from tree analysis results

import pandas as pd
import sys
from pathlib import Path

def create_combined_summary(tree_analysis_dir):
    print("=== Creating Combined Summary ===")
    
    tree_analysis_path = Path(tree_analysis_dir)
    
    # Find all summary CSV files
    summary_files = list(tree_analysis_path.rglob("*summary_*.csv"))
    
    print(f"Found {len(summary_files)} summary files")
    
    if not summary_files:
        print("No summary files found")
        return
    
    all_data = []
    
    for file_path in summary_files:
        try:
            # Read CSV
            data = pd.read_csv(file_path)
            
            if data.empty:
                continue
                
            # Add metadata
            config_name = file_path.parent.name
            data['config'] = config_name
            data['source_file'] = str(file_path)
            
            all_data.append(data)
            print(f"  Added {len(data)} rows from {file_path.name}")
            
        except Exception as e:
            print(f"  Error reading {file_path.name}: {e}")
    
    if not all_data:
        print("No valid data to combine")
        return
    
    # Combine all data
    combined_data = pd.concat(all_data, ignore_index=True)
    
    # Save combined summary
    output_file = tree_analysis_path / "all_results_summary.csv"
    combined_data.to_csv(output_file, index=False)
    
    print(f"✓ Combined summary saved: {output_file}")
    print(f"Total rows: {len(combined_data)}")
    print(f"Configurations: {combined_data['config'].nunique()}")
    print(f"Templates: {combined_data['template_name'].nunique()}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_combined_summary.py <tree_analysis_dir>")
        sys.exit(1)
    
    create_combined_summary(sys.argv[1])