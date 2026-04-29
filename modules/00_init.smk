# ============================================================================
# MODULE: 00_init.smk
# Purpose: Global variables, helper functions, and configuration
# ============================================================================

import os
from glob import glob
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import re

# Load configuration
configfile: "config.yaml"

# Global variables
GROUPS = ["target", "base"]
GENOME = ["target_base"]

def get_queries():
    """Load all query names from the config."""
    if "queries" in config:
        return list(config["queries"].keys())
    else:
        print("WARNING: No 'queries' defined in config. Using empty list.")
        return []

QUERIES = get_queries()

def get_domains(query):
    """Get domains for a specific query."""
    from glob import glob
    import os
    pattern = config["output_domains"].format(query=query)
    return [os.path.splitext(os.path.basename(f))[0] for f in glob(pattern)]

def get_pair_names(query):
    """Get pair names for a specific query."""
    import os
    pairs_file = config["filtered_pairs"].format(query=query)
    if os.path.exists(pairs_file):
        pairs = []
        with open(pairs_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    pairs.append(f"{parts[0]}_vs_{parts[1]}")
        return pairs
    else:
        print(f"Warning: {pairs_file} not found, returning empty list.")
        return []
