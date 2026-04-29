# ============================================================================
# MODULE: 12_crispr.smk
# Purpose: CRISPR guide RNA discovery and scoring
# ============================================================================

rule discover_and_score_grna:
    input:
        mrna = config["mrna"],
        idx = config["flash_fry_index"],
        header = config["header"]
    output:
        discovery = config["gRNA_discovery"],
        scored = config["crispr"],
    params:
        db_prefix = config["flash_fry_index"]
    shell:
        """
        # Step A: Discover candidates in the mRNA sequences
        java -Xmx4g -jar preset/FlashFry.jar discover \
            --database {params.db_prefix} \
            --fasta {input.mrna} \
            --output {output.discovery}

        # Step B: Score candidates (Doench for efficiency, Hsu for off-targets)
        # 'dangerous' flag checks for poly-T (terminators) and high GC content
        java -Xmx4g -jar preset/FlashFry.jar score \
            --input {output.discovery} \
            --output {output.scored} \
            --database {params.db_prefix} \
            --scoringMetrics doench2014ontarget,hsu2013,dangerous
        """
