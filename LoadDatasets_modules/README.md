# LoadDatasets_modules Documentation

## Overview

`LoadDatasets_modules/` contains 10 modularized Snakemake rules that orchestrate data acquisition, preprocessing, and indexing for the reactive bioinformatics pipeline.

**Architecture**: Modular design separates concerns into focused, reusable rule sets. The main `LoadDatasets.smk` includes all modules via `include` statements.

**Total Lines**: ~950 lines across 10 modules (vs. ~700 lines monolithic original)

## Module Structure

### Phase 1: Initialization & Configuration
**Module: `00_init_data.smk`** (~40 lines)
- Global configuration and helper functions
- Defines `groups` (target, base) and `GENOME` for analysis pairs
- Provides `get_queries()`, `get_domains()`, `get_pair_names()` helper functions

### Phase 2: Data Download
**Module: `01_download.smk`** (~160 lines)
- **Rule: `download_and_extract_data`**
  - Downloads genome, proteome, and annotation from NCBI datasets API
  - Includes robust retry logic (5 attempts) with 10-second delays
  - Performs GFF3 processing:
    - Sorts by chromosome/start/type priority/length
    - Collapses multi-part (trans-spliced) gene entries
    - Extracts mRNA sequences from CDS annotations
  - **Inputs**: taxon_id from config
  - **Outputs**: genome FASTA, proteins FASTA, annotation GFF3, mRNA FASTA

### Phase 3: GFF3 Processing
**Module: `02_gff_processing.smk`** (~90 lines)
- **Rule: `filter_gff`** - AWK: Keep only gene/mRNA/CDS features
- **Rule: `gff_to_bed`** - Convert GFF3 to BED format
- **Rule: `parse_gff`** - Extract CDS coordinates for MCScanX
  - Collapses multi-CDS per gene into min/max coordinates
  - Output: TSV with chromosome, gene_id, start, end
- **Rule: `merge_gff`** - Combine target + base GFF files
  - Creates unified GFF for synteny analysis

### Phase 4: Annotation Indexing
**Module: `03_annotation_db.smk`** (~50 lines)
- **Rule: `build_gff_db`** - Create gffutils SQLite database
  - Enables fast feature queries without re-parsing GFF
  - Parameters: merge_strategy="merge", keep_order=True
- **Rule: `samtools_index`** - Create FAI index
  - Provides random-access to genome sequences

### Phase 5: Sequence Extraction
**Module: `04_cds_extraction.smk`** (~80 lines)
- **Rule: `extract_all_cds`** (target genome)
  - Parses GFF3 CDS features
  - Reconstructs full sequences from exon fragments
  - Handles reverse-complement for "-" strand
- **Rule: `extract_all_base_cds`** (base/reference genome)
  - Same logic as above for reference species

### Phase 6: Proteome Preparation
**Module: `05_proteome_prep.smk`** (~70 lines)
- **Rule: `make_proteome_db`** - Create BLAST protein database
  - Uses `makeblastdb` for NCBI BLAST format
- **Rule: `clean_fasta`** - Reformat proteome for MCScanX
  - Maps protein IDs to gene IDs using GFF3
  - Creates cleaned FASTA with MCScanX-compatible headers

### Phase 7: Alignment & BLAST
**Module: `06_alignment.smk`** (~60 lines)
- **Rule: `diamond_make_db`** - Build DIAMOND BLAST database
- **Rule: `diamond_blastp`** - Run DIAMOND protein search
  - Parameters: e-value=1e-3, max-target-seqs=500
  - Mode: --more-sensitive for better accuracy
  - **Threads**: 8 (configurable)
- **Rule: `combine_blast`** - Merge bidirectional BLAST results
  - Combines target_vs_base + base_vs_target

### Phase 8: Synteny Analysis
**Module: `07_synteny.smk`** (~130 lines)
- **Rule: `mcscanx`** - Identify syntenic blocks
  - Parameters: s=2, g=500, m=10, e=0.01
  - **Threads**: 8
  - Output: Collinearity file with block structure
- **Rule: `duplicate_gene_classifier`** - Classify duplicated genes
  - Uses MCScanX collinearity output
- **Rule: `extract_all_pairs`** - Parse gene pairs from collinearity
  - Output: TSV with gene1, gene2 columns
  - Validates that parser found alignment blocks
- **Rule: `extract_wgd_pairs`** - Extract WGD-event duplicates
  - Groups pairs by duplication block ID

### Phase 9: Ortholog Detection (NEW)
**Module: `08_orthofinder.smk`** (~100 lines)
- **Rule: `orthofinder`** - Detect orthologous proteins
  - Tool: OrthoFinder with DIAMOND sequence search
  - **Inputs**: Target and base proteome FASTA files
  - **Outputs**: Orthogroups.tsv and other ortholog assignments
  - **Parameters**:
    - method: "diamond" (faster) or "blast" (more compatible)
    - threads: 8 (configurable)
  - **Algorithm**: Bidirectional BLAST + MCL clustering

- **Rule: `extract_orthogroups`** - Filter to 1-to-1 orthologs
  - Reads Orthogroups.tsv from OrthoFinder
  - Keeps only single-copy genes present exactly once in each genome
  - **Output**: TSV with target_gene, base_gene, orthogroup_id columns
  - Used downstream for comparative genomics analysis

### Phase 10: Auxiliary Indexing
**Module: `09_indexing.smk`** (~60 lines)
- **Rule: `gmap_build_db`** - Build GMAP genomic index
  - Enables fast mRNA-to-genome mapping
  - Used in MainPipeline for annotation extraction
- **Rule: `flashfry_index`** - Build FlashFry CRISPR index
  - Parameters: enzyme=spcas9ngg (SpCas9 NGG PAM)
  - Enables off-target search for guide RNAs

## Dependency Flow Diagram

```
01_download
    ├─> 02_gff_processing
    │       ├─> 03_annotation_db
    │       └─> 04_cds_extraction
    ├─> 05_proteome_prep
    │       └─> 06_alignment
    │           └─> 07_synteny
    ├─> 08_orthofinder (NEW)
    └─> 09_indexing
```

## Configuration Keys Used

```yaml
# From config.yaml - expected keys:
genome_group: "data/genome/{group}/genome.fna"
proteome_group: "data/proteome/{group}/proteins.faa"
annotation_group: "data/genome/{group}/annotation.gff3"
mrna_group: "data/genome/{group}/mrna.fasta"

filtered_group: "data/genome/{group}/annotation.filtered.gff3"
annotation_db_group: "data/genome/{group}/annotation.db"

target_cds: "data/genome/target/cds.fasta"
base_cds: "data/genome/base/cds.fasta"

clean_faa: "data/mcscanx/{group}.clean.faa"
protein_db_group: "data/proteome/{group}/proteins.done"

gmap_marker: "data/gmap_db/target/db.done"
flash_fry_index: "data/flashfry_index/genome"

combo_gff: "data/mcscanx/target_base.gff"
combo_collin: "data/mcscanx/target_base.collinearity"
group_pairs: "data/mcscanx/{genome}.wgd_pairs.tsv"
```

## Usage

### Run entire data loading pipeline:
```bash
snakemake -s LoadDatasets.smk -j 8
```

### Run specific module:
```bash
# Download only
snakemake -s LoadDatasets.smk rule download_and_extract_data -j 2

# GFF processing
snakemake -s LoadDatasets.smk rule parse_gff -j 1

# BLAST + synteny
snakemake -s LoadDatasets.smk rule diamond_blastp rule mcscanx -j 8

# OrthoFinder
snakemake -s LoadDatasets.smk rule orthofinder rule extract_orthogroups -j 8
```

### Dry-run to validate DAG:
```bash
snakemake -s LoadDatasets.smk -n --dag | dot -Tpng > dag.png
```

## Debugging Tips

### If download fails:
- Check internet connection
- Verify taxon_ids in config.yaml (target=4530, base=3702 for reference)
- OrthoFinder failures: Ensure both proteome directories exist with *.faa files
- Check `data/blast/OrthoFinder/Results/` for raw output

### If GFF processing fails:
- Verify GFF3 format (9 columns, tab-separated)
- Check for invalid chromosome names (spaces, special characters)
- Some GFF3 files have carriage returns; use `dos2unix` if needed

### If MCScanX fails:
- Requires both GFF and BLAST in same prefix: `data/mcscanx/target_base.*`
- Check that GFF IDs match BLAST subject/query names exactly
- MCScanX creates multiple output files; check permissions in data/mcscanx/

### Missing tool errors:
- `datasets`: Install NCBI datasets (`conda install ncbi-datasets-cli`)
- `gmap_build`: Install GMAP (`conda install gmap`)
- `orthofinder`: Install OrthoFinder (`conda install orthofinder`)
- `diamond`: Install DIAMOND (`conda install diamond`)

## Adding New Modules

To add a new module (e.g., `10_custom.smk`):

1. Create `LoadDatasets_modules/10_custom.smk` with your rules
2. Add `include: "LoadDatasets_modules/10_custom.smk"` to `LoadDatasets.smk`
3. Add output targets to `rule all` in `LoadDatasets.smk`
4. Update this README with module description
5. Test with dry-run: `snakemake -s LoadDatasets.smk -n`

## Integration with MainPipeline

After LoadDatasets.smk completes:

1. **OrthoFinder results** (`data/blast/orthofinder_1to1_pairs.tsv`) can be used in MainPipeline for:
   - Comparative genomics analysis
   - Gene family tracking across genomes
   - Ortholog-based functional annotation transfer

2. **MCScanX output** provides synteny information for:
   - Visualization of conserved gene order
   - Identification of chromosomal rearrangements
   - Whole-genome duplication (WGD) analysis

3. **Extracted sequences and indexes** enable:
   - Fast GMAP mRNA mapping in annotation module
   - CRISPR guide RNA design via FlashFry
   - Phylogenetic and domain analysis

## Performance Characteristics

| Module | Rule | Runtime (typical) | Threads | Memory |
|--------|------|-------------------|---------|--------|
| 01 | download_and_extract_data | 10-30 min | 1 | 2 GB |
| 02 | filter_gff | <1 min | 1 | <1 GB |
| 03 | build_gff_db | 2-5 min | 1 | 1-2 GB |
| 06 | diamond_blastp | 10-20 min | 8 | 4-8 GB |
| 07 | mcscanx | 5-15 min | 8 | 2-4 GB |
| 08 | orthofinder | 30-60 min | 8 | 8-16 GB |
| 09 | flashfry_index | 10-30 min | 1 | 4 GB |

**Total typical runtime**: 60-180 minutes for full pipeline (mostly download + indexing)

## Version History

- **v2.0** (Current): Modularized from monolithic (~700 → ~950 lines)
  - Added OrthoFinder module for ortholog detection
  - Improved error handling and logging
  - Better code organization and maintainability

- **v1.0**: Original monolithic LoadDatasets.smk (~700 lines)
