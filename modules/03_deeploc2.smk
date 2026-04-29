# ============================================================================
# MODULE: 03_deeploc2.smk
# Purpose: Subcellular localization prediction with DeepLoc2
# ============================================================================

rule run_deeploc2:
    input:
        fasta=config["secondary_query"]
    output:
        results=config["deeploc2_results"]
    shell:
        """
        mkdir -p $(dirname {output.results})
        deeploc2 -f {input.fasta} -o $(dirname {output.results})
        # Rename the timestamped output to the expected filename
        mv $(dirname {output.results})/results_*.csv {output.results}
        """
