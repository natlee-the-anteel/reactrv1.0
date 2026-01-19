import os
from glob import glob
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import re
configfile: "config.yaml"
GROUPS = ["target", "base"]
GENOME= ["target_base"]

def get_domains():
    from glob import glob
    import os
    return [os.path.splitext(os.path.basename(f))[0] for f in glob(config["output_domains"])]
domains = get_domains()

def get_pair_names():
    import os
    pairs_file = config["filtered_pairs"]
    if os.path.exists(pairs_file):
        pairs = []
        with open(pairs_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    pairs.append(f"{parts[0]}_vs_{parts[1]}")
        return pairs
    else:
        print(f"Warning: {pairs_file} not found, returning empty list.")
        return []
pairs = get_pair_names()

rule all:
    input:
        config["setup_marker"],
        config["secondary_query"],
        config["base_blast"],
        config["target_blast"],
        config["hmmer_scan"],
        config["hmmer_done"],
        config["sort_marker"],
        config["sort_text"],
        config["unique_pairs"],
        config["true_query"],
        
        # Dynamic per-domain sorted outputs
        expand(config["meme_xml"], domain = domains),
        expand(config["msa_aligned"], domain = domains),
        expand(config["tree"], domain = domains),
        expand(config["mrna"], domain=domains),
        expand(config["first_property"], domain=domains),
        expand(config["gmap_annotation"], domain=domains),
        
        expand(config["primers"], domain=domains, ptype = ["cloning", "validation", "qpcr"]),

        expand(config["annotations"], domain = domains),
        
        expand(config["select_genes"], domain=domains),
        expand(config["promoters"], domain = domains),
        expand(config["promoter_meme"], domain = domains),
        
        #expand(config["primersearch_txt"], domain= domains, ptype = ["cloning", "validation", "qpcr"]),
        #expand(config["primersearch_sim"], domain=domains, ptype=["cloning", "validation", "qpcr"]),
        expand(config["crispr"], domain=domains),
        
        config["pairs_to_run"],
        expand(config["pair_cds"], pair=pairs),
        expand(config["cds_align"], pair=pairs),
        expand(config["ka/ks"], pair=pairs),
        expand(config["chromosome_map"], domain=domains),
        expand(config["gene_structure"], domain=domains),
        expand(config["extra_properties"], domain=domains),

rule generate_query_fasta:
    output:
        # This reads the target filename from your config
        config["true_query"]
    run:
        # Open the output file in write mode ('w')
        with open(output[0], "w") as f:
            # Write the string stored in the config
            f.write(config["query_contents"])

rule setup_output:
    output:
        touch(config["setup_marker"])
    shell:
        "mkdir -p output && touch {output}"
rule extract_filtered_target_hits: 
    input:
        query=config["true_query"],
        db=config["target_proteome"],
    output:
        filtered=config["secondary_query"]
    params:
        db_prefix=config["target_prot_prefix"],
        tmp_blast=config["target_tmp"],
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

        # -----------------------------
        # Load sequences
        # -----------------------------
        query_seqs = SeqIO.to_dict(SeqIO.parse(input.query_fasta, "fasta"))
        query_ids = set(query_seqs)

        combined_seqs = {
            rec.id: rec for rec in SeqIO.parse(input.combined_fasta, "fasta")
        }

        # -----------------------------
        # Parse HMMER tblout
        # domain → gene → best isoform
        # -----------------------------
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

        # -----------------------------
        # Clean old outputs
        # -----------------------------
        for f in glob(os.path.join(domain_dir, "*.fasta")):
            os.remove(f)

        # -----------------------------
        # Write per-domain FASTAs
        # -----------------------------
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

rule meme_per_domain:
    input:
        fasta = config["domain_fasta"]
    output:
        xml = config["meme_xml"]
    params:
        nmotifs = 10,
        minw = 6,
        maxw = 50,
        meme_opts = "-protein -mod zoops"
    shell:
        """
        # 1. Get the directory path dynamically from the output file
        #    (e.g., if output is "results/meme/dom1/meme.xml", this becomes "results/meme/dom1")
        OUT_DIR=$(dirname {output.xml})

        # 2. Create the directory (safely handles if it doesn't exist)
        mkdir -p "$OUT_DIR"

        # 3. Run MEME using the dynamic directory variable
        meme {input.fasta} {params.meme_opts} \
             -nmotifs {params.nmotifs} \
             -minw {params.minw} \
             -maxw {params.maxw} \
             -oc "$OUT_DIR"
        """

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
rule extract_mrna_for_domain:
    input:
        domain_fasta=config["domain_fasta"],
        mrna_fasta=config["target_mrna"],
        gff_db=config["target_annotation"]
    output:
        mrna=config["mrna"]
    run:
        import os
        import gffutils
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.mrna), exist_ok=True)
        db = gffutils.FeatureDB(input.gff_db)

        # 1. Get Protein IDs from HMMER (e.g., XXG55436.1)
        # We need these to match the headers in mrna.fasta
        target_ids = {rec.id.split()[0] for rec in SeqIO.parse(input.domain_fasta, "fasta")}
        
        print(f"[INFO] {wildcards.domain}: Searching for {len(target_ids)} target IDs.")

        # 2. Extract mRNA sequences from FASTA
        # Based on your grep, the mrna.fasta uses XXG IDs directly.
        # We don't even need the GFF bridge if the Protein IDs match the mRNA FASTA IDs!
        mrna_records = []
        for rec in SeqIO.parse(input.mrna_fasta, "fasta"):
            clean_id = rec.id.split()[0]
            if clean_id in target_ids:
                mrna_records.append(rec)
        
        # 3. Validation and Writing
        if not mrna_records:
            print(f"[WARNING] No matches found for {wildcards.domain}!")
            # Diagnostic: Print one ID from each to see why they don't match
            if target_ids:
                sample_target = list(target_ids)[0]
                # Re-parse one record to get a sample
                sample_fasta = next(SeqIO.parse(input.mrna_fasta, "fasta")).id
                print(f"[DEBUG] Sample Target ID: {sample_target}")
                print(f"[DEBUG] Sample mRNA FASTA ID: {sample_fasta}")
        
        SeqIO.write(mrna_records, output.mrna, "fasta")
        print(f"[SUCCESS] Wrote {len(mrna_records)} sequences to {output.mrna}")
rule design_primers_per_domain:
    input:
        mrna_fasta=lambda wildcards: config["mrna"].format(domain=wildcards.domain)
    output:
        primer_file=config["primers"]
    wildcard_constraints:
        ptype="|".join(["cloning", "validation", "qpcr"])
    params:
        product_ranges={
            "qpcr":"80-150",
            "validation":"300-800",
            "cloning":"700-1500"
        }
    run:
        import os
        from Bio import SeqIO
        from subprocess import run, PIPE

        os.makedirs(os.path.dirname(output.primer_file), exist_ok=True)
        product_range = params.product_ranges[wildcards.ptype]

        seqs = list(SeqIO.parse(input.mrna_fasta, "fasta"))
        if not seqs:
            print(f"Warning: no mRNA sequences for domain {wildcards.domain} → primers will not be designed")
            open(output.primer_file, "w").close()
            return

        # Clear file
        open(output.primer_file, 'w').close()

        for record in seqs:
            seq_id = record.id
            sequence = str(record.seq).replace("U", "T")

            primer_input = "\n".join([
                f"SEQUENCE_ID={seq_id}",
                f"SEQUENCE_TEMPLATE={sequence}",
                "PRIMER_TASK=generic",
                "PRIMER_OPT_SIZE=20",
                "PRIMER_MIN_SIZE=18",
                "PRIMER_MAX_SIZE=25",
                "PRIMER_MIN_TM=57.0",
                "PRIMER_OPT_TM=60.0",
                "PRIMER_MAX_TM=63.0",
                f"PRIMER_PRODUCT_SIZE_RANGE={product_range}",
                "P3_FILE_FLAG=0",
                "PRIMER_EXPLAIN_FLAG=1",
                "="
            ]) + "\n"

            result = run("primer3_core", input=primer_input.encode(), stdout=PIPE, stderr=PIPE)
            with open(output.primer_file, "a") as out_f:
                out_f.write(f"== {seq_id} ==\n")
                out_f.write(result.stdout.decode() + "\n")
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

rule compute_protein_properties_per_domain:
    input:
        fasta=config["domain_fasta"]
    output:
        props=config["first_property"]
    run:
        import os
        from Bio import SeqIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        os.makedirs(os.path.dirname(output.props), exist_ok=True)

        with open(output.props, "w") as out_f:
            out_f.write("ID\tMW\tpI\n")

            for rec in SeqIO.parse(input.fasta, "fasta"):
                seq = str(rec.seq).replace("*", "").replace("X", "")
                if not seq or len(seq) < 2:
                    continue
                try:
                    analysis = ProteinAnalysis(seq)
                    mw = analysis.molecular_weight()
                    pi = analysis.isoelectric_point()
                    out_f.write(f"{rec.id}\t{mw:.2f}\t{pi:.2f}\n")
                except Exception:
                    out_f.write(f"{rec.id}\tNA\tNA\n")

rule extract_annotations_by_domain:
    input:
        gmap_gff=config["gmap_annotation"],
        gff_db=config["target_annotation"]
    output:
        gff=config["annotations"]
    run:
        import gffutils
        import os
        import sys
        import re

        target_domain = wildcards.domain
        target_protein_ids = set()

        sys.stderr.write(f"\n[STEP 1] Extracting IDs from GMAP file for: {target_domain}\n")

        # 1. Parse the GMAP GFF3 for Protein IDs
        # We look for 'Name=XP_XXXXXX.1' or 'Target=XP_XXXXXX.1'
        if os.path.exists(input.gmap_gff):
            with open(input.gmap_gff, 'r') as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    # Use regex to find XP/NP IDs (e.g., XP_015624862.1)
                    matches = re.findall(r'[XN]P_\d+\.\d+', line)
                    for m in matches:
                        target_protein_ids.add(m)
        
        sys.stderr.write(f"[STEP 1] Found {len(target_protein_ids)} IDs: {list(target_protein_ids)[:3]}...\n")

        if not target_protein_ids:
            with open(output.gff, 'w') as f: f.write("##gff-version 3\n")
            return

        # 2. Database Lookup
        db = gffutils.FeatureDB(input.gff_db)
        relevant_feature_ids = set()
        
        # Iterate CDS features to find matching protein_ids
        for cds in db.features_of_type("CDS"):
            pids_in_gff = cds.attributes.get('protein_id', [])
            if any(p in target_protein_ids for p in pids_in_gff):
                relevant_feature_ids.add(cds.id)
                
                # Get the Family Tree (Parent mRNA, Grandparent Gene)
                for parent in db.parents(cds):
                    relevant_feature_ids.add(parent.id)
                    for gp in db.parents(parent):
                        relevant_feature_ids.add(gp.id)
                
                # Get Siblings (Exons)
                for parent in db.parents(cds):
                    for child in db.children(parent):
                        relevant_feature_ids.add(child.id)

        sys.stderr.write(f"[STEP 2] Matched {len(relevant_feature_ids)} features in Database.\n")

        # 3. Write Output
        os.makedirs(os.path.dirname(output.gff), exist_ok=True)
        with open(output.gff, "w") as out:
            out.write("##gff-version 3\n")
            if relevant_feature_ids:
                # Sort by chromosome and start position
                sorted_ids = sorted(relevant_feature_ids, 
                                  key=lambda x: (db[x].chrom, db[x].start))
                for fid in sorted_ids:
                    out.write(str(db[fid]) + "\n")

        sys.stderr.write(f"[STEP 3] Success! Wrote {output.gff}\n\n")
rule extract_selected_genes:
    input:
        gmap_gff=config["gmap_annotation"]
    output:
        txt=config["select_genes"]
    shell:
        # This looks for "Name=" in the GFF, pulls the ID, and removes duplicates
        """
        grep "Name=" {input.gmap_gff} | sed 's/.*Name=//;s/;.*//' | sort -u > {output.txt}
        """

rule extract_promoters:
    input:
        genome=config["target_genome"],
        gff_db=config["target_annotation"], # Use the DB we already built
        genes=config["select_genes"]
    output:
        promoters=config["promoters"]
    run:
        from Bio import SeqIO
        import gffutils
        import sys

        # 1. Load Genome
        sys.stderr.write("Loading genome into memory...\n")
        genome_dict = SeqIO.to_dict(SeqIO.parse(input.genome, "fasta"))
        
        # 2. Connect to GFF database
        db = gffutils.FeatureDB(input.gff_db)
        
        # 3. Read target IDs
        with open(input.genes) as f:
            target_ids = {line.strip() for line in f if line.strip()}
        
        promoter_len = 2000
        output_records = []

        sys.stderr.write(f"Processing {len(target_ids)} genes...\n")

        # 4. Extract sequences
        for pid in target_ids:
            # Find the CDS feature that has this protein_id
            found_feature = None
            for cds in db.features_of_type("CDS"):
                if pid in cds.attributes.get('protein_id', []):
                    # We want the FIRST CDS (start codon)
                    # For + strand, it's the one with the smallest 'start'
                    # For - strand, it's the one with the largest 'end'
                    if not found_feature:
                        found_feature = cds
                    else:
                        if cds.strand == '+' and cds.start < found_feature.start:
                            found_feature = cds
                        elif cds.strand == '-' and cds.end > found_feature.end:
                            found_feature = cds
            
            if not found_feature:
                sys.stderr.write(f"Warning: Could not find CDS for {pid}\n")
                continue

            chrom = found_feature.chrom
            strand = found_feature.strand
            
            if strand == '+':
                p_start = max(1, found_feature.start - promoter_len)
                p_end = found_feature.start - 1
            else:
                p_start = found_feature.end + 1
                p_end = min(len(genome_dict[chrom]), found_feature.end + promoter_len)

            # Slice and handle reverse complement for minus strand
            promo_seq = genome_dict[chrom].seq[p_start-1:p_end]
            if strand == '-':
                promo_seq = promo_seq.reverse_complement()

            header = f">{pid}_{chrom}_{p_start}_{p_end}_{strand}"
            output_records.append(f"{header}\n{promo_seq}\n")

        # 5. Write final file
        with open(output.promoters, "w") as f:
            f.writelines(output_records)
            
        sys.stderr.write(f"Done! Wrote {len(output_records)} promoters to {output.promoters}\n")
rule meme_promoters:
    input:
        fasta=config["promoters"]
    output:
        xml=config["promoter_xml"],
        html=config["promoter_meme"],
        skipped=config["promoter_skip"] # Keep this
    shell:
        """
        # Count sequences in the fasta
        seq_count=$(grep -c ">" {input.fasta} || true)

        if [ "$seq_count" -lt 1 ]; then
            echo "Too few sequences ($seq_count). Skipping MEME."
            mkdir -p $(dirname {output.xml})
            touch {output.skipped}
            # Create dummy files so Snakemake is happy
            touch {output.xml} {output.html}
        else
            meme {input.fasta} -dna -oc $(dirname {output.xml}) -mod zoops -nmotifs 10 -minw 6 -maxw 50 -revcomp
            # CRITICAL: Create the SKIPPED file as a placeholder even if we didn't skip
            touch {output.skipped}
        fi
        """

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

rule primersearch:
    input:
        genome = config["target_genome"],
        primers = config["primersearch_txt"]
    output:
        report = config["primersearch_sim"]
    params:
        mismatch = 0,       # max mismatches allowed
    shell:
        """
        primersearch -seqall {input.genome} \
                    -infile {input.primers} \
                    -mismatch {params.mismatch} \
                    -outfile {output.report}
        """
rule primers_to_primersearch:
    input:
        primer3_txt=config["primers"]
    output:
        primersearch=config["primersearch_txt"]
    params:
        valid_types=["cloning", "validation", "qpcr"]
    run:
        import re

        if wildcards.ptype not in params.valid_types:
            raise ValueError(f"Invalid primer type: {wildcards.ptype}")

        out_lines = []

        with open(input.primer3_txt) as f:
            content = f.read()

        # Split by primer blocks labeled by "== geneid =="
        blocks = re.split(r"==\s*(\S+)\s*==", content)

        for i in range(1, len(blocks), 2):
            seq_id = blocks[i].strip()
            block = blocks[i + 1]

            # Find all PRIMER_LEFT_n and PRIMER_RIGHT_n entries
            lefts = re.findall(r"PRIMER_LEFT_(\d+)_SEQUENCE=([A-Za-z]+)", block)
            rights = re.findall(r"PRIMER_RIGHT_(\d+)_SEQUENCE=([A-Za-z]+)", block)

            # Pair by shared index
            pairs = {}
            for idx, seq in lefts:
                pairs[idx] = {"left": seq}
            for idx, seq in rights:
                if idx in pairs:
                    pairs[idx]["right"] = seq

            # Write output in EMBOSS primersearch format
            for idx in sorted(pairs.keys(), key=int):
                pair = pairs[idx]
                if "left" in pair and "right" in pair:
                    line = f"{seq_id}_pair{idx} {pair['left']} {pair['right']}"
                    out_lines.append(line)

        with open(output.primersearch, "w") as out_f:
            out_f.write("\n".join(out_lines) + "\n")

rule calculate_conservation_and_kaks_input:
    input:
        wgd=config["mcscanx_wgd_pairs"],
        unique=config["unique_pairs"]
    output:
        filtered=config["filtered_pairs"],
        stats=config["conservation_stats"],
        pair_list=config["pairs_to_run"]
    run:
        import os
        
        def clean(s):
            # Removes versions (.1), colons (0:), and whitespace
            return s.strip().split('.')[0].rstrip(':').upper()

        # 1. Build a "Family Registry" from your BLAST results
        # We store every gene that BLAST says belongs to your family
        family_registry = set()
        with open(input.unique) as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    family_registry.add(clean(parts[0]))
                    family_registry.add(clean(parts[1]))
        
        print(f"DEBUG: Family Registry contains {len(family_registry)} unique Gene IDs.")

        # 2. Filter MCScanX: Keep pair if BOTH genes are in the Registry
        # This confirms that the syntenic block contains your actual genes of interest
        conserved_pairs = []
        with open(input.wgd) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2: continue
                
                g1_orig, g2_orig = parts[0], parts[1]
                g1_clean, g2_clean = clean(g1_orig), clean(g2_orig)
                
                if g1_clean in family_registry and g2_clean in family_registry:
                    conserved_pairs.append((g1_orig, g2_orig))

        # 3. Write outputs
        os.makedirs(os.path.dirname(output.filtered), exist_ok=True)
        with open(output.filtered, "w") as out_tsv, open(output.pair_list, "w") as out_list:
            out_tsv.write("gene1\tgene2\n")
            for g1, g2 in conserved_pairs:
                out_tsv.write(f"{g1}\t{g2}\n")
                # Format for Snakemake wildcards {pair}
                out_list.write(f"{g1}_vs_{g2}\n")

        # 4. Global Stats for the paper
        total_pairs = len(conserved_pairs)
        with open(output.stats, "w") as f:
            f.write(f"Species-Agnostic Conservation Analysis\n")
            f.write(f"Total Family Members in Registry: {len(family_registry)}\n")
            f.write(f"Pairs maintained in Syntenic Blocks: {total_pairs}\n")
        
        print(f"SUCCESS: Found {total_pairs} conserved pairs.")
rule extract_pair_cds:
    input:
        cds_base = config["base_cds"],
        cds_target = config["target_cds"]
    output:
        config["pair_cds"]
    run:
        from Bio import SeqIO  # <--- import inside run
        import pandas as pd

        # Load all CDS sequences into a dictionary
        all_cds = {}
        for f in [input.cds_base, input.cds_target]:
            all_cds.update(SeqIO.to_dict(SeqIO.parse(f, "fasta")))

        # Split wildcard to get gene IDs
        gene1, gene2 = wildcards.pair.split("_vs_")
        records = []
        for g in [gene1, gene2]:
            if g in all_cds:
                records.append(all_cds[g])
            else:
                print(f"Warning: {g} not found in CDS files")

        if records:
            SeqIO.write(records, output[0], "fasta")
rule codon_align_pair:
    input:
        cds_fasta = config["pair_cds"]
    output:
        prot_fasta = config["pair_cds_fasta"],
        prot_aln = config["cds_aligned_fasta"],
        codon_aln = config["cds_align"]
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import tempfile
        import os

        # Create output directory
        os.makedirs(os.path.dirname(output.prot_fasta), exist_ok=True)

        # Translate CDS to protein (stop at first stop codon)
        records = []
        for rec in SeqIO.parse(input.cds_fasta, "fasta"):
            prot_rec = SeqRecord(rec.seq.translate(to_stop=True), id=rec.id, description="")
            records.append(prot_rec)

        # Write temporary protein FASTA for alignment
        tmp_prot = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta")
        SeqIO.write(records, tmp_prot.name, "fasta")
        tmp_prot.close()

        # Save protein FASTA as final output
        SeqIO.write(records, output.prot_fasta, "fasta")

        # Run MUSCLE 5 (pairwise alignment)
        shell(f"muscle -align {tmp_prot.name} -output {output.prot_aln} > /dev/null 2>&1")

        # Check alignment produced
        if not os.path.exists(output.prot_aln) or os.path.getsize(output.prot_aln) == 0:
            raise RuntimeError(f"MUSCLE failed to produce alignment for {input.cds_fasta}")

        # Codon alignment with Pal2Nal
        shell(f"pal2nal.pl {output.prot_aln} {input.cds_fasta} -output fasta > {output.codon_aln}")
rule codon_to_axt:
    input:
        codon = config["cds_align"]
    output:
        axt = config["axt_convert"]
    run:
        from Bio import SeqIO
        records = list(SeqIO.parse(input.codon, "fasta"))
        if len(records) != 2:
            raise ValueError(f"{input.codon} does not have exactly 2 sequences")
        with open(output.axt, "w") as f:
            f.write(f"1 {records[0].id} {records[1].id}\n")
            f.write(str(records[0].seq) + "\n")
            f.write(str(records[1].seq) + "\n")
rule calculate_kaks:
    input:
        axt = config["axt_convert"]
    output:
        kaks = config["ka/ks"]
    shell:
        """
        mkdir -p $(dirname {output.kaks})
        KaKs_Calculator -i {input.axt} -o {output.kaks} -m NG
        """

rule discover_and_score_grna:
    input:
        mrna = config["mrna"],
        idx = config["flash_fry_index"],
        header = config["header"]
    output:
        discovery = config["gRNA_discovery"],
        scored = config["crispr"],
    params:
        db_prefix = config["flash_fry_index"]
    shell:
        """
        # Step A: Discover candidates in the mRNA sequences
        java -Xmx4g -jar preset/FlashFry.jar discover \
            --database {params.db_prefix} \
            --fasta {input.mrna} \
            --output {output.discovery}

        # Step B: Score candidates (Doench for efficiency, Hsu for off-targets)
        # 'dangerous' flag checks for poly-T (terminators) and high GC content
        java -Xmx4g -jar preset/FlashFry.jar score \
            --input {output.discovery} \
            --output {output.scored} \
            --database {params.db_prefix} \
            --scoringMetrics doench2014ontarget,hsu2013,dangerous
        """

rule plot_chromosomes:
    input:
        gff = config["target_gff"],
        candidates = config["domain_fasta"],
        fai = config["target_fai"]
    output:
        config["chromosome_map"]
    run:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from pygenomeviz import GenomeViz
        from BCBio import GFF

        # --- Helper: Recursive Search ---
        def feature_matches(feature, targets):
            attributes_to_check = [
                feature.id,
                feature.qualifiers.get("Name", [""])[0],
                feature.qualifiers.get("protein_id", [""])[0],
                feature.qualifiers.get("locus_tag", [""])[0]
            ]
            for attr in attributes_to_check:
                for target in targets:
                    if target and attr and target in attr:
                        return True, target
            for sub in feature.sub_features:
                match, matched_id = feature_matches(sub, targets)
                if match:
                    return True, matched_id
            return False, None

        # 1. Load Candidates for this domain
        target_genes = set()
        with open(input.candidates) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    target_genes.add(line[1:].split()[0])
                elif not any(c in line for c in "ATCGN") and line: 
                    target_genes.add(line)

        print(f"DEBUG: Processing domain {wildcards.domain} - {len(target_genes)} IDs.")

        # 2. Read Chromosome Lengths
        fai_df = pd.read_csv(input.fai, sep="\t", header=None, names=["name", "length", "x", "y", "z"])
        seqid2size = dict(zip(fai_df["name"], fai_df["length"]))

        # 3. Setup Plot
        gv = GenomeViz()
        gv.set_scale_bar()

        # 4. Parse GFF
        found_features = []
        with open(input.gff) as in_handle:
            for rec in GFF.parse(in_handle):
                if rec.id in seqid2size:
                    for feature in rec.features:
                        is_match, matched_id = feature_matches(feature, target_genes)
                        if is_match:
                            found_features.append({
                                "chrom": rec.id,
                                "start": int(feature.location.start),
                                "end": int(feature.location.end),
                                "strand": feature.strand if feature.strand is not None else 1,
                                "id": matched_id
                            })

        # 5. Add Tracks
        chroms_with_hits = set(f["chrom"] for f in found_features)
        
        if not chroms_with_hits:
            print(f"WARNING: No matches found for domain {wildcards.domain}.")
            # Plot empty map with first chromosome as fallback
            first_chrom = list(seqid2size.keys())[0]
            gv.add_feature_track(first_chrom, seqid2size[first_chrom])
        else:
            for seqid, size in seqid2size.items():
                if seqid in chroms_with_hits:
                    track = gv.add_feature_track(seqid, size)
                    for f in found_features:
                        if f["chrom"] == seqid:
                            track.add_feature(f["start"], f["end"], f["strand"], 
                                              label=f["id"], plotstyle="arrow", facecolor="red")

        # 6. Save
        gv.savefig(output[0])
rule plot_gene_structure:
    input:
        gff = config["target_gff"],
        candidates = config["domain_fasta"]
    output:
        directory(config["gene_structure"])
    run:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from BCBio import GFF
        from dna_features_viewer import GraphicFeature, GraphicRecord
        import os

        # 1. Create Output Directory
        os.makedirs(output[0], exist_ok=True)

        # 2. Load Candidates for this domain
        target_genes = set()
        with open(input.candidates) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    target_genes.add(line[1:].split()[0])
                elif not any(c in line for c in "ATCGN") and line:
                    target_genes.add(line)

        # 3. Helper: Recursive Search
        def feature_matches(feature, targets):
            attributes_to_check = [
                feature.id,
                feature.qualifiers.get("Name", [""])[0],
                feature.qualifiers.get("protein_id", [""])[0],
                feature.qualifiers.get("locus_tag", [""])[0]
            ]
            for attr in attributes_to_check:
                for target in targets:
                    if target and attr and target in attr:
                        return True, target
            for sub in feature.sub_features:
                match, matched_id = feature_matches(sub, targets)
                if match:
                    return True, matched_id
            return False, None

        # 4. Helper: Extract Subfeatures
        def get_all_subfeatures(feature):
            subs = []
            if feature.type == "CDS":
                subs.append(feature)
            for sub in feature.sub_features:
                subs.extend(get_all_subfeatures(sub))
            return subs

        # 5. Parse GFF and Plot
        print(f"DEBUG: Generating plots for domain {wildcards.domain} ({len(target_genes)} targets)...")
        
        with open(input.gff) as in_handle:
            for rec in GFF.parse(in_handle):
                for feature in rec.features:
                    is_match, matched_id = feature_matches(feature, target_genes)
                    
                    if is_match:
                        # Found a target
                        print(f"DEBUG: Plotting {matched_id}...")
                        
                        parts = get_all_subfeatures(feature)
                        if not parts:
                            parts = [feature]

                        graphic_features = []
                        gene_start = int(feature.location.start)
                        gene_end = int(feature.location.end)
                        gene_len = gene_end - gene_start

                        for part in parts:
                            start_rel = int(part.location.start) - gene_start
                            end_rel = int(part.location.end) - gene_start
                            strand = part.strand if part.strand is not None else 1
                            
                            graphic_features.append(GraphicFeature(
                                start=start_rel, end=end_rel, strand=strand,
                                color="#ffcccc", label="CDS" if part.type == "CDS" else part.type
                            ))

                        record = GraphicRecord(sequence_length=gene_len, features=graphic_features)
                        ax, _ = record.plot(figure_width=10)
                        ax.set_title(f"Gene Structure: {matched_id}", loc='left')
                        
                        safe_name = matched_id.replace(".", "_")
                        outfile = os.path.join(output[0], f"{safe_name}.png")
                        plt.savefig(outfile, bbox_inches='tight')
                        plt.close()

rule calculate_extra_properties:
    input:
        config["domain_fasta"]
    output:
        config["extra_properties"]
    run:
        from Bio import SeqIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        import pandas as pd

        data = []

        # Parse the fasta file
        for record in SeqIO.parse(input[0], "fasta"):
            # Clean sequence: Remove stop (*) and handle ambiguity (X)
            # ProtParam strictly requires standard Amino Acids. 
            # We strip non-standard chars for calculation purposes.
            seq_str = str(record.seq).replace("*", "")
            clean_seq = "".join([aa for aa in seq_str if aa in "ACDEFGHIKLMNPQRSTVWY"])

            if not clean_seq:
                continue

            try:
                # Initialize ProtParam analysis
                analyser = ProteinAnalysis(clean_seq)
                
                # 1. Molecular Weight (in Daltons -> convert to kDa)
                mw = analyser.molecular_weight()
                
                # 2. Isoelectric Point (pI)
                pi = analyser.isoelectric_point()
                
                # 3. Instability Index (<40 stable, >40 unstable)
                instability = analyser.instability_index()
                
                # 4. GRAVY (Grand Average of Hydropathy)
                gravy = analyser.gravy()

                # 5. Aliphatic Index 
                # Formula: X(Ala) + 2.9*X(Val) + 3.9*(X(Ile) + X(Leu)) 
                # (using mole percentages)
                aa_percents = analyser.get_amino_acids_percent() # returns fraction (e.g. 0.10)
                aliphatic = (aa_percents.get('A', 0) + 
                             2.9 * aa_percents.get('V', 0) + 
                             3.9 * (aa_percents.get('I', 0) + aa_percents.get('L', 0))) * 100

                # 6. Amino Acid Count (Total length)
                length = len(clean_seq)

                data.append({
                    "GeneID": record.id,
                    "Length": length,
                    "MW_kDa": round(mw / 1000, 2),
                    "pI": round(pi, 2),
                    "Instability_Index": round(instability, 2),
                    "Aliphatic_Index": round(aliphatic, 2),
                    "GRAVY": round(gravy, 3),
                    "Stability_Pred": "Stable" if instability < 40 else "Unstable"
                })

            except Exception as e:
                print(f"WARNING: Could not process {record.id}: {e}")

        # Save to CSV
        if data:
            df = pd.DataFrame(data)
            df.to_csv(output[0], index=False)
        else:
            # Create empty file if no data to prevent Snakemake errors
            with open(output[0], 'w') as f:
                f.write("GeneID,Length,MW_kDa,pI,Instability_Index,Aliphatic_Index,GRAVY,Stability_Pred\n")

