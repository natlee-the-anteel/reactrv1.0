# ============================================================================
# MODULE: 04_hmmer.smk
# Purpose: Domain detection using HMMER and Pfam database
# ============================================================================

rule run_hmmer_on_query:
    input:
        fasta=config["secondary_query"],
        pfam=config["pfam"]
    output:
        tbl=config["sort_text"],
        domtbl=config["domain_tbl"],
        hits=config["hmmer_hits"],
        skipped=config["hmmer_skip"]
    threads: 4
    run:
        import os
        import subprocess
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.tbl), exist_ok=True)

        records = list(SeqIO.parse(input.fasta, "fasta"))

        def skip(reason):
            for f in [output.tbl, output.domtbl, output.hits]:
                with open(f, "w") as out:
                    out.write(f"# HMMER skipped: {reason}\n")
            with open(output.skipped, "w") as out:
                out.write(reason + "\n")
            print(f"[HMMER] Skipping: {reason}")

        # --- HARD GUARDS ---
        if not records:
            skip("query.fasta is empty")
            return

        valid = [
            r for r in records
            if len(str(r.seq).replace("*", "").replace("X", "")) >= 20
        ]

        if not valid:
            skip("no valid protein sequences after filtering")
            return

        # Rewrite cleaned FASTA (important)
        SeqIO.write(valid, input.fasta, "fasta")

        subprocess.run(
            [
                "hmmscan",
                "--cpu", str(threads),
                "--tblout", output.tbl,
                "--domtblout", output.domtbl,
                input.pfam,
                input.fasta,
            ],
            check=True
        )

        # Extract hit sequences (safe even if no hits)
        hit_ids = set()
        with open(output.domtbl) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                hit_ids.add(line.split()[0])

        hits = [r for r in valid if r.id in hit_ids]

        SeqIO.write(hits, output.hits, "fasta")

        if not hits:
            skip("no Pfam hits detected")

rule combine_hits_and_query:
    input:
        hits=config["hit_seqs"],
        query=config["secondary_query"]
    output:
        combined=config["combo_hits"]
    run:
        from Bio import SeqIO
        with open(output.combined, "w") as out_f:
            for f in [input.query, input.hits]:
                for rec in SeqIO.parse(f, "fasta"):
                    SeqIO.write(rec, out_f, "fasta")

rule run_hmmer:
    input:
        fasta=config["combo_hits"],
        pfam_db=config["pfam"]
    output:
        out=config["hmmer_scan"]
    shell:
        """
        hmmscan --cpu 8 --tblout {output.out} {input.pfam_db} {input.fasta}
        """

rule split_by_domain:
    input:
        hmmer_tbl=config["hmmer_scan"],
        combined_fasta=config["combo_hits"],
        query_fasta=config["secondary_query"]
    output:
        done=config["hmmer_done"]
    params:
        max_total=8,
        max_bg=6,
    run:
        import os
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from glob import glob

        def gene_key(seq_id):
            base = seq_id.split(".")[0]
            return re.sub(r"(X\d+)$", "", base)

        domain_dir = os.path.dirname(output.done)
        os.makedirs(domain_dir, exist_ok=True)

        # Load sequences
        query_seqs = SeqIO.to_dict(SeqIO.parse(input.query_fasta, "fasta"))
        query_ids = set(query_seqs)

        combined_seqs = {
            rec.id: rec for rec in SeqIO.parse(input.combined_fasta, "fasta")
        }

        # Parse HMMER tblout
        # domain → gene → best isoform
        domain_hits = defaultdict(dict)

        with open(input.hmmer_tbl) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.split()
                domain = cols[0]
                seq_id = cols[2]
                bitscore = float(cols[5])

                gkey = gene_key(seq_id)

                if (
                    gkey not in domain_hits[domain]
                    or bitscore > domain_hits[domain][gkey][1]
                ):
                    domain_hits[domain][gkey] = (seq_id, bitscore)

        # Clean old outputs
        for f in glob(os.path.join(domain_dir, "*.fasta")):
            os.remove(f)

        # Write per-domain FASTAs
        for domain, gene_map in domain_hits.items():

            hits = sorted(
                gene_map.values(),
                key=lambda x: x[1],
                reverse=True
            )

            query_hits = [sid for sid, _ in hits if sid in query_ids]
            bg_hits = [sid for sid, _ in hits if sid not in query_ids]

            if not query_hits:
                continue

            bg_hits = bg_hits[:params.max_bg]
            final_ids = (query_hits + bg_hits)[:params.max_total]

            if len(final_ids) < 2:
                continue

            domain_file = os.path.join(domain_dir, f"{domain}.fasta")
            with open(domain_file, "w") as out_f:
                for sid in final_ids:
                    SeqIO.write(combined_seqs[sid], out_f, "fasta")

        with open(output.done, "w") as f:
            f.write("done\n")

rule split_query_by_domain:
    input:
        hmmer_tbl=config["sort_text"],
        combined_fasta=config["hmmer_hits"]
    output:
        done=config["sort_marker"]
    run:
        import os
        from collections import defaultdict
        from Bio import SeqIO

        domain_dir = os.path.dirname(output.done)
        os.makedirs(domain_dir, exist_ok=True)

        # Parse HMMER table
        domain_to_seqids = defaultdict(set)
        with open(input.hmmer_tbl) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                domain = parts[0]
                seq_id = parts[2]
                domain_to_seqids[domain].add(seq_id)

        # Read all sequences
        seq_records = SeqIO.to_dict(SeqIO.parse(input.combined_fasta, "fasta"))

        # Clean up any old FASTAs
        for fname in os.listdir(domain_dir):
            if fname.endswith(".fasta"):
                os.remove(os.path.join(domain_dir, fname))

        # Write only domains with 2+ sequences
        for domain, seq_ids in domain_to_seqids.items():
            valid_ids = [seq_id for seq_id in seq_ids if seq_id in seq_records]
            if len(valid_ids) < 2:
                continue
            with open(os.path.join(domain_dir, f"{domain}.fasta"), "w") as out_f:
                for seq_id in valid_ids:
                    SeqIO.write(seq_records[seq_id], out_f, "fasta")

        with open(output.done, "w") as f:
            f.write("done\n")
