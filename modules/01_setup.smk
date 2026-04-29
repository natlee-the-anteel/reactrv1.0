# ============================================================================
# MODULE: 01_setup.smk
# Purpose: Setup and initialization rules
# ============================================================================

rule setup_output:
    output:
        touch(config["setup_marker"])
    shell:
        "rm -rf outputs/ && mkdir -p $(dirname {output}) && touch {output}"

rule generate_query_fasta:
    output:
        config["true_query"]
    run:
        with open(output[0], "w") as f:
            f.write(config["queries"][wildcards.query])
