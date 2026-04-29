# ============================================================================
# MAIN PIPELINE (Modularized)
# Purpose: Include all modular workflows and define final aggregation rule
# ============================================================================

# Include initialization module (must be first)
include: "modules/00_init.smk"

# Include all functional modules
include: "modules/01_setup.smk"
include: "modules/02_query_filtering.smk"
include: "modules/03_deeploc2.smk"
include: "modules/04_hmmer.smk"
include: "modules/05_alignment_phylo.smk"
include: "modules/06_domain_analysis.smk"
include: "modules/07_annotation.smk"
include: "modules/08_gmap.smk"
include: "modules/09_primers.smk"
include: "modules/10_conservation.smk"
include: "modules/11_visualization.smk"
include: "modules/12_crispr.smk"

# ============================================================================
# MAIN AGGREGATION RULE
# ============================================================================

rule all:
    input:
        # Setup and query generation (per query)
        expand(config["setup_marker"], query=QUERIES),
        expand(config["secondary_query"], query=QUERIES),
        expand(config["true_query"], query=QUERIES),
        expand(config["deeploc2_results"], query=QUERIES),
        
        # BLAST outputs (per query)
        expand(config["base_blast"], query=QUERIES),
        expand(config["target_blast"], query=QUERIES),
        expand(config["unique_pairs"], query=QUERIES),
        
        # HMMER outputs (per query)
        expand(config["hmmer_scan"], query=QUERIES),
        expand(config["hmmer_done"], query=QUERIES),
        expand(config["sort_marker"], query=QUERIES),
        expand(config["sort_text"], query=QUERIES),
        
        # Per-domain outputs: Compute dynamically per query
        *[expand(config["meme_xml"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["msa_aligned"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["tree"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["mrna"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["first_property"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["gmap_annotation"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["primers"], query=q, domain=get_domains(q), ptype=["cloning", "validation", "qpcr"]) for q in QUERIES],
        *[expand(config["annotations"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["select_genes"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["promoters"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["promoter_meme"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["crispr"], query=q, domain=get_domains(q)) for q in QUERIES],
        
        # Pairs/conservation analysis
        *[expand(config["pairs_to_run"], query=q) for q in QUERIES],
        *[expand(config["pair_cds"], query=q, pair=get_pair_names(q)) for q in QUERIES],
        *[expand(config["cds_align"], query=q, pair=get_pair_names(q)) for q in QUERIES],
        *[expand(config["ka/ks"], query=q, pair=get_pair_names(q)) for q in QUERIES],
        *[expand(config["chromosome_map"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["gene_structure"], query=q, domain=get_domains(q)) for q in QUERIES],
        *[expand(config["extra_properties"], query=q, domain=get_domains(q)) for q in QUERIES],
