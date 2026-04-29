# LoadDatasets_modules/05_proteome_prep.smk
# Purpose: Prepare proteome sequences and databases for analysis

rule make_proteome_db:
    """Create BLAST protein database from proteome FASTA.
    
    Builds NCBI BLAST-format database for quick similarity searches.
    """
    input: config["proteome_group"]
    output: config["protein_db_group"]
    threads: 4
    shell:
        """
        makeblastdb -in {input} -dbtype prot -out data/proteome/{wildcards.group}/proteins
        touch {output}
        """

rule clean_fasta:
    """Clean and reformat proteome FASTA for MCScanX analysis.
    
    Maps protein IDs to gene IDs using GFF3 annotations.
    Ensures headers are compatible with MCScanX format requirements.
    """
    input:
        fasta=lambda wc: f"data/proteome/{wc.group}/proteins.faa",
        gff=lambda wc: f"data/genome/{wc.group}/annotation.filtered.gff3"
    output:
        clean="data/mcscanx/{group}.clean.faa"
    run:
        import os
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.clean), exist_ok=True)

        # Step 1: Build protein → gene mapping from GFF
        prot_to_gene = {}
        with open(input.gff) as gff_file:
            for line in gff_file:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attrs = dict(x.split("=", 1) for x in parts[8].split(";") if "=" in x)
                gene_id = attrs.get("ID") or attrs.get("Name") or attrs.get("gene")
                parent  = attrs.get("Parent")
                if gene_id and parent:
                    prot_to_gene[parent] = gene_id

        print(f"[{wildcards.group}] Found {len(prot_to_gene)} protein→gene mappings")

        # Step 2: Rewrite FASTA headers
        written = 0
        with open(output.clean, "w") as out_f:
            for rec in SeqIO.parse(input.fasta, "fasta"):
                new_id = prot_to_gene.get(rec.id, rec.id)
                rec.id = new_id
                rec.description = ""
                SeqIO.write(rec, out_f, "fasta")
                written += 1

        if written == 0:
            raise ValueError(
                f"[{wildcards.group}] clean_fasta produced ZERO sequences — "
                "check FASTA/GFF ID compatibility"
            )

        print(f"[{wildcards.group}] Wrote {written} sequences to {output.clean}")
