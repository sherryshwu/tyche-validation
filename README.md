# Computational validation pipeline for TyCHE
A computational pipeline for evaluating the performance of TyCHE (Type-linked Clocks for Heterogeneous Evolution) that infers time-resolved phylogeneies of distinctively evolving cell populations and benchmarking tree topology and ancestral cell type inference against the ground truth using simulated data generated with SimBLE.

## Overview
We perform simulation-based studies of B cell evolution by:
1) Simulating B-cell receptor (BCR) data with various evolutionary scenarios
2) Running Bayesian phylogenetic analysis using BEAST 2 and Dowser
3) Analyzing reconstruction accuracy of trees and ancestral cell types
4) Generate publication-ready figures and summaries

## Project structure
```
scripts/
├── 00_master_pipeline.sh
├── 01_simulate_data.sh
├── 02_submit_beast_phases.sh
├── 03_tree_analysis_main.sh
├── analysis
│   ├── create_combined_summary.py
│   ├── tree_analysis.R
│   └── tree_functions.R
├── beast
│   ├── create_beast_job_combinations.sh
│   ├── generate_configs.py
│   ├── run_beast_dowser.R
│   └── run_beast_dowser.sh
└── plotting
    ├── publication_plots_main.R
    └── tree_plotting.R
```

## Usage
### Complete pipeline
Running the complete analysis pipeline:
`bash scripts/00_master_pipeline.sh [simulation_name]`

## Directory Structure
Files are organized in:
[simulation_name]/
├── data/                     
│   ├── raw/
│   ├── processed/
├── results/
│   ├── main_analysis/
│   │   ├── irrev/
│   │   └── rev/
│   ├── differentiation_analysis/  
│   │   ├── irrev/
│   ├── sub_analysis/
│   │   ├── irrev/                 
├── logs/
└── configs/