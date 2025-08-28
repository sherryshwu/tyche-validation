#!/usr/bin/env python3
import json
import sys
import os

def generate_config_files(simulation_dir):
    # Create config directory within the simulation
    config_dir = os.path.join(simulation_dir, "configs")
    os.makedirs(config_dir, exist_ok=True)
    
    print(f"Creating config files in: {config_dir}")
    
    # Base configuration template
    base_config = {
        "LOCATIONS": [
            {
                "name": {"__enum__": "LocationName.GC"},
                "sample_times": list(range(10, 201, 10)),
                "mutation_rate": 1.0,
                "max_population": 1000,
                "migration_rate": 2,
                "sample_size": 10
            },
            {
                "name": {"__enum__": "LocationName.OTHER"},
                "sample_times": list(range(10, 201, 10)),
                "mutation_rate": 0.0,
                "max_population": 1000,
                "migration_rate": 0,
                "sample_size": 10
            }
        ],
        "HEAVY_SHM_PER_SITE": 0.0008908272571108565,
        "LIGHT_SHM_PER_SITE": 0.0004923076923076923,
        "TARGET_MUTATIONS_HEAVY": 10,
        "TARGET_MUTATIONS_LIGHT": 0,
        "UNIFORM": False,
        "SELECTION": True,
        "RESULTS_DIR": "/home/sample_config_results",
        "MULTIPLIER": 2,
        "DEV": False,
        "FASTA": False,
        "VERBOSE": False,
        "CDR_DIST": "exponential",
        "CDR_VAR": 0.995,
        "FWR_DIST": "constant",
        "FWR_VAR": 2
    }
    
    # Create descriptive filename with all parameters
    def create_filename(ratio_label, evolution_type):
        return f"config_ratio_{ratio_label}_{evolution_type}.json"
    
    config_count = 0
    
    # Define sample times variations
    sample_times_options = [
        (list(range(50, 201, 50)), "50to200")
    ]
    
    # Vary sample size ratios combined with selection and sample times
    sample_ratios = [
        (12, 12, "1to1"),  # 1:1 ratio
        (6, 18, "1to3"),   # 1:3 ratio
    ]
    
    # Define evolution types with their corresponding UNIFORM and SELECTION settings
    evolution_types = [
        {
            "type": "sel",
            "uniform": False,
            "selection": True,
            "description": "Selective evolution (UNIFORM=False, SELECTION=True)"
        },
        {
            "type": "neu",
            "uniform": True,
            "selection": False,
            "description": "Neutral evolution (UNIFORM=True, SELECTION=False)"
        }
    ]
    
    for sample_times, sample_times_label in sample_times_options:
        for gc_size, other_size, ratio_label in sample_ratios:
            for evolution in evolution_types:
                config = json.loads(json.dumps(base_config))
                
                # Update sample times for both locations
                config["LOCATIONS"][0]["sample_times"] = sample_times
                config["LOCATIONS"][1]["sample_times"] = sample_times
                
                # Update sample sizes
                config["LOCATIONS"][0]["sample_size"] = gc_size
                config["LOCATIONS"][1]["sample_size"] = other_size
                
                # Update evolution parameters
                config["UNIFORM"] = evolution["uniform"]
                config["SELECTION"] = evolution["selection"]
                
                filename = create_filename(sample_times_label, ratio_label, evolution["type"])
                
                with open(os.path.join(config_dir, filename), 'w') as f:
                    json.dump(config, f, indent=2)
                
                config_count += 1
                print(f"Generated: {filename}")
                print(f"  {evolution['description']}")
    
    print(f"\nGenerated {config_count} config files in {config_dir}")
    print(f" - Sample ratios x Evolution types x Sample times: {config_count} configs")
    print(f"   ({len(sample_ratios)} ratios x {len(evolution_types)} evolution types x {len(sample_times_options)} sample times)")
    print(f" - Total: {config_count} configs")
    
    # Print summary of evolution types
    print(f"\nEvolution types generated:")
    for evolution in evolution_types:
        print(f" - {evolution['type']}: {evolution['description']}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_configs.py <simulation_dir>")
        sys.exit(1)
    
    simulation_dir = sys.argv[1]
    
    print(f"Arguments received:")
    print(f"  simulation_dir: {simulation_dir}")
    
    generate_config_files(simulation_dir)