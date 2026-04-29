# Snakemake Pipeline Modules

This directory contains modularized Snakemake rules organized by functional domains. Each module focuses on a specific aspect of the bioinformatics analysis pipeline.

## Module Overview

### **00_init.smk** - Initialization & Configuration
- Global variables and helper functions
- Query and domain name resolution
- Configuration loading
- **Must be included first**

### **01_setup.smk** - Setup & Initialization
- Output directory initialization
- Cleanup of previous analysis outputs
- Query FASTA generation
- Setup markers for dependency tracking

### **02_query_filtering.smk** - Query Filtering & BLAST
- BLAST against target and base genomes
- Hit filtering by e-value and identity
- Unique pair identification
- Hit sequence extraction

### **03_deeploc2.smk** - Subcellular Localization
- DeepLoc2 subcellular localization prediction
- CSV output generation

### **04_hmmer.smk** - Domain Detection (HMMER/Pfam)
- HMMER scan against Pfam database
- Hit combination and merging
- Domain-specific FASTA splitting
- Query-specific domain extraction

### **05_alignment_phylo.smk** - Alignment & Phylogenetics
- Multiple sequence alignment (MUSCLE)
- Phylogenetic tree inference (IQ-TREE)
- Tree quality validation

### **06_domain_analysis.smk** - Domain-Specific Analysis
- Motif discovery (MEME)
- Protein property calculation (MW, pI, instability, GRAVY)
- Physicochemical properties

### **07_annotation.smk** - Genomic Annotation & Gene Selection
- mRNA extraction for domains
- Gene annotation extraction
- Promoter sequence extraction (2kb upstream)
- Promoter motif analysis (MEME)

### **08_gmap.smk** - GMAP Genomic Mapping
- mRNA to genome mapping (GFF3 output)
- Gene-to-coordinate resolution

### **09_primers.smk** - Primer Design & Search
- Primer3 primer design (cloning, validation, qPCR)
- EMBOSS PrimerSearch for primer-target validation
- Primer evaluation

### **10_conservation.smk** - Conservation & Evolution
- Gene pair identification from MCScanX results
- Codon-aware sequence alignment (Pal2Nal)
- Ka/Ks calculation (non-synonymous/synonymous substitution rates)
- Conservation statistics

### **11_visualization.smk** - Data Visualization
- Chromosome-level gene localization maps
- Gene structure diagrams (exon/CDS visualization)
- Publication-quality figures

### **12_crispr.smk** - CRISPR Guide RNA Design
- Guide RNA candidate discovery (FlashFry)
- Off-target and efficiency scoring
- Ranked gRNA tables

## Module Dependency Flow

```
00_init
  ↓
01_setup
  ├→ 02_query_filtering
  │   └→ 03_deeploc2
  │
  ├→ 04_hmmer → 05_alignment_phylo
  │                └→ 06_domain_analysis
  │
  ├→ 07_annotation → 08_gmap
  │                    └→ 09_primers
  │                        └→ 12_crispr
  │
  └→ 10_conservation
```

## Usage

The main pipeline is invoked from `MainPipeline.smk`, which includes all modules:

```bash
snakemake -s MainPipeline.smk [options]
```

To run a specific analysis module:

```bash
snakemake -s MainPipeline.smk rule_name [options]
```

## Best Practices

1. **Module Independence**: Each module is self-contained with complete documentation
2. **Clear Dependencies**: Use config variables to link modules
3. **Error Handling**: Modules include graceful handling of missing/empty data
4. **Logging**: Diagnostic messages for troubleshooting
5. **Scalability**: Wildcards support multiple queries and domains

## Adding New Modules

To add a new analysis module:

1. Create `XX_name.smk` in the `modules/` directory
2. Document module purpose and dependencies
3. Add include statement to `MainPipeline.smk` in logical order
4. Add outputs to the `rule all` input list
5. Add configuration entries to `config.yaml`

## Debugging

To debug a specific module:

```bash
snakemake -s MainPipeline.smk rule_name --dry-run --verbose
```

To inspect module dependencies:

```bash
snakemake -s MainPipeline.smk --dag | dot -Tpdf > dag.pdf
```
