#!/usr/bin/env snakemake
# LoadDatasets.smk
# Purpose: Modularized data loading and preprocessing pipeline
# 
# This pipeline orchestrates:
# 1. Downloading reference genomes and annotations from NCBI
# 2. Processing and filtering GFF3 annotations
# 3. Building indexes for downstream analysis
# 4. Preparing proteomes and extracting sequences
# 5. Running BLAST and MCScanX for synteny analysis
# 6. Detecting orthologs with OrthoFinder
#
# Module structure: 10 focused modules in LoadDatasets_modules/

import os
from subprocess import run

configfile: "config.yaml"

# Include all modularized rules
include: "LoadDatasets_modules/00_init_data.smk"
include: "LoadDatasets_modules/01_download.smk"
include: "LoadDatasets_modules/02_gff_processing.smk"
include: "LoadDatasets_modules/03_annotation_db.smk"
include: "LoadDatasets_modules/04_cds_extraction.smk"
include: "LoadDatasets_modules/05_proteome_prep.smk"
include: "LoadDatasets_modules/06_alignment.smk"
include: "LoadDatasets_modules/07_synteny.smk"
include: "LoadDatasets_modules/08_orthofinder.smk"
include: "LoadDatasets_modules/09_indexing.smk"

rule all:
    """
    Final target: Ensure all data loading steps complete.
    
    Dependencies:
    1. Genomes, proteomes, annotations downloaded
    2. GFF3 processed and indexed
    3. CDS sequences extracted
    4. Protein databases created
    5. GMAP and FlashFry indexes built
    6. Synteny analysis completed (MCScanX)
    
    Note: OrthoFinder is optional (run separately if needed):
        snakemake -s LoadDatasets.smk -R orthofinder extract_orthogroups -j 8
    """
    input:
        expand(config["genome_group"], group=groups),
        expand(config["proteome_group"], group=groups),
        expand(config["annotation_group"], group=groups),
        expand(config["mrna_group"], group=groups),

        expand(config["filtered_group"], group=groups),
        expand(config["bed_group"], group=groups),

        expand(config["annotation_db_group"], group=groups),

        expand(config["protein_db_group"], group=groups),
        config["gmap_marker"],
        config["target_cds"],
        config["base_cds"],
        config["flash_fry_index"],

        expand(config["group_pairs"], genome=GENOME),
        config["combo_gff"],
        config["genome_fna_fai"]
