# LoadDatasets_modules/06_alignment.smk
# Purpose: Perform bidirectional BLAST alignments for ortholog identification

rule diamond_make_db:
    """Create DIAMOND BLAST database from cleaned proteome FASTA."""
    input:
        fasta=config["clean_faa"]
    output:
        db=config["blast_dmnd"]
    threads: 4
    shell:
        """
        mkdir -p data/blast
        if [ ! -f {output.db} ]; then
            diamond makedb --in {input.fasta} -d {output.db}
        fi
        """

rule diamond_blastp:
    """Run DIAMOND BLAST for fast protein similarity search.
    
    Compares query proteome against reference database.
    Parameters:
    - e-value: 1e-3 (significance threshold)
    - max-target-seqs: 500 (keep top 500 hits per query)
    - --more-sensitive: More accurate but slower than default
    """
    input:
        query=config["query_clean"],
        db=config["blast_db"]
    output:
        out=config["output_diamond"]
    threads: 8
    params:
        evalue=1e-3,
        max_target_seqs=500,
        sensitive="--more-sensitive"
    shell:
        """
        mkdir -p data/blast
        diamond blastp \
            --query {input.query} \
            --db {input.db} \
            --out {output.out} \
            --evalue {params.evalue} \
            --max-target-seqs {params.max_target_seqs} \
            --outfmt 6 \
            {params.sensitive} \
            --threads {threads}
        """

rule combine_blast:
    """Merge bidirectional BLAST results for synteny analysis.
    
    Combines:
    - target_vs_base: Target proteins queried against base database
    - base_vs_target: Base proteins queried against target database
    """
    input:
        blast1=config["target_vs_base"],
        blast2=config["base_vs_target"]
    output:
        combined=config["combo_blast"]
    shell:
        """
        mkdir -p data/mcscanx
        cat {input.blast1} {input.blast2} > {output.combined}
        """
