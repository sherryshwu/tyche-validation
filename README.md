# Computational validation pipeline for TyCHE
A computational pipeline for evaluating the performance of TyCHE (Type-linked Clocks for Heterogeneous Evolution) that infers time-resolved phylogeneies of distinctively evolving cell populations and benchmarking tree topology and ancestral cell type inference against the ground truth using simulated data generated with SimBLE.

## Overview
We perform simulation-based studies of B cell evolution by:
1) Simulating B-cell receptor (BCR) data with various evolutionary scenarios
2) Running Bayesian phylogenetic analysis using BEAST 2 and Dowser
3) Analyzing reconstruction accuracy of trees and ancestral cell types
4) Generate publication-ready figures and summaries

## Repository structure
```
scripts
├── analysis
│   ├── calculate_tree_metrics.R
│   ├── consolidate_tree_metrics.py
│   └── summarize_convergence.R
├── phylogenetics
│   ├── create_phylo_job_grid.sh
│   ├── generate_simble_configs.py
│   ├── run_beast_dowser.R
│   └── submit_beast_dowser.sh
├── pipeline
│   ├── 00_run_pipeline.sh
│   ├── 01_simulate_data.sh
│   ├── 02_submit_beast_phases.sh
│   ├── 03_analyze_trees.sh
│   └── 04_generate_figures.sh
├── setup
│   └── import_gc_reentry_data.sh
├── utils
│   └── phylo_utilities.R
└── visualization
    ├── create_differentiation_figures.R
    ├── create_main_figures.R
    └── plot_trees.R
```

## Usage
### Complete pipeline
Running the complete analysis pipeline:
```
bash scripts/00_run_pipeline.sh [simulation_name]
```

### Individual steps
1) Data simulation: `sbatch scripts/01_simulate_data.sh [simulation_name]`

2) Beast running: `bash scripts/02_submit_beast_phases.sh [simulation_name] [analysis_type] [rev_suffix]`

3) Tree analysis: `bash scripts/03_analyze_trees.sh [simulation_name] [analysis_type] [rev_suffix]`

4) Figures generation: `bash scripts/04_generate_figures.sh [analysis_type]`

## Analysis Types
The pipeline supports the following analysis types:
- main_analysis: core phylogenetic analysis across all evolutionary configurations
- differentiation_analysis: applying advanced three-state phylogenetic model to reconstruct the timing of cell differentiation
- sub_analysis*: subgroup analysis performed on various sampled timepoints to investigate the phenomenon of time-dependent rate decay

## File Organization
Files are organized in:
```
[simulation_name]/
├── data/                     
│   ├── raw/
│   ├── processed/
├── figures/
│   ├── main_analysis/
│   ├── differentiation_analysis/
├── results/
│   ├── main_analysis/
│   │   ├── irrev/
│   │   └── rev/
│   ├── differentiation_analysis/  
│   │   └── irrev/             
├── logs/
└── configs/
```