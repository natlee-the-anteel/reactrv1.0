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
        nseq = len(records)

        unique_seqs = {
            str(r.seq).replace("-", "")
            for r in records
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
                    "-s", input.aln,
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
