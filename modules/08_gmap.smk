# ============================================================================
# MODULE: 08_gmap.smk
# Purpose: GMAP genomic mapping of mRNA to genes
# ============================================================================

rule gmap_per_domain:
    input:
        mrna=config["mrna"]
    output:
        gff=config["gmap_annotation"]
    params:
        db_dir=config["target_gmap"],
        db_name="db"
    shell:
        """
        mkdir -p $(dirname {output.gff})

        # Check if the input file has data (-s means file exists and size > 0)
        if [ -s {input.mrna} ]; then
            # We call 'gmap' directly now, relying on the Conda environment
            gmap -D {params.db_dir} -d {params.db_name} \
                 -f gff3_match_cdna -n 0 -t 4 {input.mrna} > {output.gff}
        else
            # Create a dummy/empty GFF if there are no inputs to prevent errors
            echo "## No mRNA sequences to map for domain" > {output.gff}
        fi
        """
