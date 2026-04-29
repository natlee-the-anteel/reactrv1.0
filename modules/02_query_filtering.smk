# ============================================================================
# MODULE: 02_query_filtering.smk
# Purpose: Query filtering and BLAST analysis
# ============================================================================

rule extract_filtered_target_hits: 
    input:
        query=config["true_query"],
        db=config["target_proteome"],
    output:
        filtered=config["secondary_query"]
    params:
        db_prefix=config["target_prot_prefix"],
        tmp_blast=lambda wc: config["target_tmp"].format(query=wc.query),
        evalue_thresh=1e-5,
        min_identity=30,   # percent identity threshold
        max_hits=5,           # max hits per Arabidopsis query
    run:
        import os
        from Bio import SeqIO
        from collections import defaultdict

        os.makedirs(os.path.dirname(params.tmp_blast), exist_ok=True)

        # Build BLAST DB if needed
        if not os.path.exists(f"{params.db_prefix}.pin"):
            shell(f"makeblastdb -in {input.db} -dbtype prot -out {params.db_prefix}")

        # Run BLASTP
        shell(f"blastp -query {input.query} -db {params.db_prefix} "
              f"-out {params.tmp_blast} -outfmt 6 -evalue {params.evalue_thresh} -max_target_seqs 1000")

        # Parse and filter
        hits_by_query = defaultdict(list)
        with open(params.tmp_blast) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue
                qid, sid = cols[0], cols[1]
                identity = float(cols[2])
                evalue = float(cols[10])
                if evalue <= float(params.evalue_thresh) and identity >= float(params.min_identity):
                    hits_by_query[qid].append((sid, evalue, identity))

        # Limit to top N hits per query
        selected_ids = set()
        for qid, hit_list in hits_by_query.items():
            sorted_hits = sorted(hit_list, key=lambda x: x[1])  # sort by evalue
            for hit in sorted_hits[:params.max_hits]:
                selected_ids.add(hit[0])

        # Extract protein sequences
        records = []
        for rec in SeqIO.parse(input.db, "fasta"):
            if rec.id in selected_ids:
                records.append(rec)

        os.makedirs(os.path.dirname(output.filtered), exist_ok=True)
        with open(output.filtered, "w") as out_f:
            SeqIO.write(records, out_f, "fasta")

rule blastp_base:
    input:
        query=config["secondary_query"],
        db=config["base_proteome"],
    output:
        out=config["base_blast"]
    params:
        db_name=config["base_prot_prefix"]
    shell:
        """
        if [ ! -f {params.db_name}.pin ]; then
            makeblastdb -in {input.db} -dbtype prot -out {params.db_name}
        fi
        blastp -query {input.query} -db {params.db_name} -out {output.out} -outfmt 6 -evalue 1e-5 -max_target_seqs 5
        """

rule convert_query_to_tsv:
    input:
        query=config["secondary_query"],
        pairs=config["filtered_pairs"]
    output:
        out=config["target_blast"]
    run:
        import os
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.out), exist_ok=True)

        with open(output.out, "w") as out_tsv:
            for record in SeqIO.parse(input.query, "fasta"):
                # Dummy TSV format matching BLAST outfmt 6:
                # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                out_tsv.write(f"{record.id}\tMockTarget\t100.00\t{len(record.seq)}\t0\t0\t1\t{len(record.seq)}\t1\t{len(record.seq)}\t0.0\t200\n")

rule extract_hit_ids:
    input:
        base=config["base_blast"],
        target=config["target_blast"]
    output:
        ids=config["hit_ids"]
    run:
        hits = set()
        for infile in [input.base, input.target]:
            with open(infile) as f:
                hits.update(line.split()[1] for line in f if line.strip())
        with open(output.ids, "w") as out:
            out.write("\n".join(sorted(hits)) + "\n")

rule extract_hit_sequences:
    input:
        ids=config["hit_ids"],
        base_fasta=config["base_proteome"],
        target_fasta=config["target_proteome"]
    output:
        fasta=config["hit_seqs"]
    run:
        from Bio import SeqIO

        # Load the list of hit IDs
        with open(input.ids) as f:
            ids = set(line.strip() for line in f)

        # Ensure output directory exists
        import os
        os.makedirs(os.path.dirname(output.fasta), exist_ok=True)

        # Open output FASTA for writing
        with open(output.fasta, "w") as out:
            for db_path in [input.base_fasta, input.target_fasta]:
                for rec in SeqIO.parse(db_path, "fasta"):
                    if rec.id in ids:
                        SeqIO.write(rec, out, "fasta")

rule extract_unique_pairs_no_dup_within_species:
    input:
        blast_tsv=config["base_blast"]
    output:
        unique_pairs=config["unique_pairs"]
    run:
        seen_target = set()
        seen_base = set()
        with open(input.blast_tsv) as infile, open(output.unique_pairs, "w") as outfile:
            for line in infile:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.strip().split('\t')
                target_gene, base_gene = cols[0], cols[1]

                if target_gene not in seen_target and base_gene not in seen_base:
                    seen_target.add(target_gene)
                    seen_base.add(base_gene)
                    outfile.write(f"{target_gene}\t{base_gene}\n")
