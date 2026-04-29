# LoadDatasets_modules/00_init_data.smk
# Purpose: Global initialization and configuration for data loading pipeline

import os
from subprocess import run

# Define analysis groups (target species vs reference/base species)
groups = ["target", "base"]
GENOME = ["target_base"]

# Helper function to list all domains from disk
def get_domains():
    """
    Scans output/{query}/domain_sorted/ to find all domain directories.
    Returns list of domain names (subdirectory names).
    Requires MainPipeline to have been run to create these directories.
    """
    domains = []
    if os.path.exists("output"):
        for item in os.listdir("output"):
            domain_dir = f"output/{item}/domain_sorted"
            if os.path.isdir(domain_dir):
                for domain in os.listdir(domain_dir):
                    domain_path = os.path.join(domain_dir, domain)
                    if os.path.isdir(domain_path) and domain != "domain_sorted":
                        if domain not in domains:
                            domains.append(domain)
    return domains if domains else ["unknown"]

# Helper function to extract all queries
def get_queries():
    """
    Returns list of all query names from config.queries.
    Used for expanding rules across multiple queries.
    """
    return list(config.get("queries", {}).keys())

# Helper function to extract pair names for WGD analysis
def get_pair_names():
    """
    Returns list of genome pair combinations for multi-genome analysis.
    Used for MCScanX and DIAMOND BLAST processing.
    """
    pairs = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i+1:]:
            pairs.append(f"{g1}_{g2}")
            pairs.append(f"{g2}_{g1}")
    return pairs
