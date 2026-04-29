# ============================================================================
# MODULE: 10_conservation.smk
# Purpose: Conservation analysis and KaKs calculations
# ============================================================================

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
        from Bio import SeqIO

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
