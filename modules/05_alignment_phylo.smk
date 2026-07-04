# ============================================================================
# MODULE: 05_alignment_phylo.smk
# Purpose: Multiple sequence alignment and phylogenetic tree inference
# ============================================================================

rule msa_per_domain:
    input:
        domain_fasta=config["domain_fasta"],
        split_done=config["hmmer_done"],
    output:
        aligned=config["msa_aligned"]
    threads: 4
    shell:
        """
        muscle -quiet -align {input.domain_fasta} -output {output.aligned}.tmp -threads {threads}

        # Replace dashes in sequence IDs with underscores to make MEME/IQ-TREE safe
        awk '/^>/ {{ gsub(/-/, "_"); }} {{ print }}' {output.aligned}.tmp > {output.aligned}

        rm {output.aligned}.tmp
        """

rule iqtree_per_domain:
    input:
        aln=config["msa_aligned"]
    output:
        tree=config["tree"],
        skipped=config["tree_skip"]
    log:
        config["tree_log"]
    threads: 4
    run:
        import os
        import subprocess
        from Bio import SeqIO

        outdir = os.path.dirname(output.tree)
        os.makedirs(outdir, exist_ok=True)

        records = list(SeqIO.parse(input.aln, "fasta"))

        def pairwise_identity(seq_a, seq_b):
            seq_a = str(seq_a)
            seq_b = str(seq_b)
            aligned = sum(1 for a, b in zip(seq_a, seq_b) if a != "-" or b != "-")
            if aligned == 0:
                return 1.0
            matches = sum(1 for a, b in zip(seq_a, seq_b) if a == b and a != "-")
            return matches / aligned

        kept_records = []
        collapsed = []
        collapse_threshold = 0.995

        for record in records:
            seq = str(record.seq)
            is_duplicate = False
            for kept in kept_records:
                if pairwise_identity(seq, str(kept.seq)) >= collapse_threshold:
                    collapsed.append((record.id, kept.id))
                    is_duplicate = True
                    break
            if not is_duplicate:
                kept_records.append(record)

        reduced_aln = os.path.join(outdir, f"{wildcards.domain}.collapsed.fasta")
        SeqIO.write(kept_records, reduced_aln, "fasta")

        nseq = len(kept_records)
        unique_seqs = {
            str(r.seq).replace("-", "")
            for r in kept_records
            if str(r.seq).replace("-", "")
        }

        def skip(reason):
            with open(output.tree, "w") as f:
                f.write(f"# IQ-TREE skipped: {reason}\n")
            with open(output.skipped, "w") as f:
                f.write(reason + "\n")
            with open(log[0], "w") as f:
                f.write(reason + "\n")
            print(f"[IQ-TREE] Skipping {wildcards.domain}: {reason}")

        # --- SKIP CONDITIONS ---
        if nseq < 2:
            skip(f"only {nseq} sequence(s)")
            return

        if len(unique_seqs) < 2:
            skip("all sequences identical or empty")
            return

        try:
            subprocess.run(
                [
                    "iqtree",
                    "-m", "TEST",
                    "-bb", "1000",
                    "-nt", str(threads),
                    "-s", reduced_aln,
                    "-pre", os.path.join(outdir, wildcards.domain),
                ],
                stdout=open(log[0], "w"),
                stderr=subprocess.STDOUT,
                check=True
            )

            # IMPORTANT: keep SKIPPED but empty it
            with open(output.skipped, "w") as f:
                f.write("")

        except subprocess.CalledProcessError:
            skip("IQ-TREE failed internally (uninformative alignment)")
