# LoadDatasets_modules/03_annotation_db.smk
# Purpose: Build annotation databases for efficient lookup

rule build_gff_db:
    """Create gffutils SQLite database from filtered GFF3.
    
    Enables fast queries of genes, mRNA, and CDS features
    without parsing the GFF file repeatedly.
    """
    input: config["filtered_group"]
    output: config["annotation_db_group"]
    threads: 1
    run:
        import gffutils
        gff_file = input[0]
        db_file = output[0]

        if os.path.exists(db_file):
            print(f"[{wildcards.group}] GFF database already exists, skipping")
            return

        print(f"[{wildcards.group}] Building GFF database from {gff_file}...")
        try:
            db = gffutils.create_db(
                gff_file,
                db_file,
                merge_strategy="merge",
                keep_order=True,
                force=True,
                verbose=False
            )
            print(f"[{wildcards.group}] Successfully built database with {len(db.conn.execute('SELECT COUNT(*) FROM features').fetchone())} features")
        except Exception as e:
            print(f"[{wildcards.group}] Error building database: {e}")
            raise

rule samtools_index:
    """Create samtools FAI index for fast random access to genome sequences."""
    input:
        config["target_genome"]
    output:
        config["genome_fna_fai"]
    shell:
        "samtools faidx {input}"
