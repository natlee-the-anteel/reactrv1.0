# LoadDatasets_modules/09_indexing.smk
# Purpose: Build auxiliary indexes for downstream analyses

rule gmap_build_db:
    """Build GMAP index for genomic mapping of mRNA sequences.
    
    GMAP (Genomic Mapping and Alignment Program) allows fast
    mapping of cDNA/EST sequences to their genomic coordinates.
    Used in MainPipeline for annotation extraction.
    """
    input: 
        genome = config["target_genome"]
    output: 
        done = config["gmap_marker"]
    threads: 4
    shell:
        """
        mkdir -p data/gmap_db/target
        # Use the conda-installed gmap_build with -L flag for large genomes
        gmap_build -L -D data/gmap_db/target -d db {input.genome}
        touch {output.done}
        """

rule flashfry_index:
    """Build FlashFry index for CRISPR guide RNA design.
    
    FlashFry is a tool for genome-wide discovery of off-target
    matches for CRISPR-Cas9 guide RNAs. This rule creates the
    indexed database needed for fast off-target searches.
    
    Parameters:
    - enzyme: spcas9ngg (SpCas9 with NGG PAM pattern)
    """
    input:
        ref = config["target_genome"]
    output:
        idx = config["flash_fry_index"],
        header = config["header"]
    shell:
        """
        mkdir -p data/flashfry_tmp
        mkdir -p data/flashfry_index
        
        java -Xmx4g -jar preset/FlashFry.jar index \
            --reference {input.ref} \
            --database data/flashfry_index/genome \
            --tmpLocation data/flashfry_tmp \
            --enzyme spcas9ngg
        """
