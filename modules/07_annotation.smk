# ============================================================================
# MODULE: 07_annotation.smk
# Purpose: Genomic annotation, mRNA extraction, and gene selection
# ============================================================================

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
        target_ids = {rec.id.split()[0] for rec in SeqIO.parse(input.domain_fasta, "fasta")}
        
        print(f"[INFO] {wildcards.domain}: Searching for {len(target_ids)} target IDs.")

        # 2. Extract mRNA sequences from FASTA
        mrna_records = []
        for rec in SeqIO.parse(input.mrna_fasta, "fasta"):
            clean_id = rec.id.split()[0]
            if clean_id in target_ids:
                mrna_records.append(rec)
        
        # 3. Validation and Writing
        if not mrna_records:
            print(f"[WARNING] No matches found for {wildcards.domain}!")
            if target_ids:
                sample_target = list(target_ids)[0]
                sample_fasta = next(SeqIO.parse(input.mrna_fasta, "fasta")).id
                print(f"[DEBUG] Sample Target ID: {sample_target}")
                print(f"[DEBUG] Sample mRNA FASTA ID: {sample_fasta}")
        
        SeqIO.write(mrna_records, output.mrna, "fasta")
        print(f"[SUCCESS] Wrote {len(mrna_records)} sequences to {output.mrna}")

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
        gff_db=config["target_annotation"],
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
        skipped=config["promoter_skip"]
    shell:
        """
        # Count sequences in the fasta
        seq_count=$(grep -c ">" {input.fasta} || true)

        if [ "$seq_count" -lt 2 ]; then
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
