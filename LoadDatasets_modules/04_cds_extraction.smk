# LoadDatasets_modules/04_cds_extraction.smk
# Purpose: Extract CDS sequences from genome and GFF annotations

rule extract_all_cds:
    """Extract all CDS sequences from target genome.
    
    Parses GFF3 CDS features, reconstructs full coding sequences
    from exon fragments, and handles reverse-complement for - strand.
    """
    input:
        genome_fasta = config["target_genome"],
        gff_file = config["filtered_target"]
    output:
        config["target_cds"]
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import re

        # Load genome
        genome = SeqIO.to_dict(SeqIO.parse(input.genome_fasta, "fasta"))

        # Parse GFF for CDS features
        cds_dict = {}
        with open(input.gff_file) as gff:
            for line in gff:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                chrom, _, _, start, end, _, strand, phase, attr = parts
                m_prot = re.search(r'protein_id=([^;]+)', attr)
                gene_id = m_prot.group(1) if m_prot else re.search(r'Parent=([^;]+)', attr).group(1).split(",")[0]
                cds_dict.setdefault(gene_id, []).append((chrom, int(start), int(end), strand))

        # Collect all CDS sequences
        all_records = []
        for gene_id, features in cds_dict.items():
            features_sorted = sorted(features, key=lambda x: x[1])
            seq_fragments = [genome[f[0]].seq[f[1]-1:f[2]] for f in features_sorted]
            full_seq = Seq("").join(seq_fragments)
            if features_sorted[0][3] == "-":
                full_seq = full_seq.reverse_complement()
            all_records.append(SeqRecord(full_seq, id=gene_id, description="CDS"))

        # Write all CDS sequences to a single FASTA
        SeqIO.write(all_records, output[0], "fasta")

rule extract_all_base_cds:
    """Extract all CDS sequences from base/reference genome."""
    input:
        genome_fasta = config["base_genome"],
        gff_file = config["filtered_based"]
    output:
        config["base_cds"]
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import re

        # Load genome
        genome = SeqIO.to_dict(SeqIO.parse(input.genome_fasta, "fasta"))

        # Parse GFF for CDS features
        cds_dict = {}
        with open(input.gff_file) as gff:
            for line in gff:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                chrom, _, _, start, end, _, strand, phase, attr = parts
                m_prot = re.search(r'protein_id=([^;]+)', attr)
                gene_id = m_prot.group(1) if m_prot else re.search(r'Parent=([^;]+)', attr).group(1).split(",")[0]
                cds_dict.setdefault(gene_id, []).append((chrom, int(start), int(end), strand))

        # Collect all CDS sequences
        all_records = []
        for gene_id, features in cds_dict.items():
            features_sorted = sorted(features, key=lambda x: x[1])
            seq_fragments = [genome[f[0]].seq[f[1]-1:f[2]] for f in features_sorted]
            full_seq = Seq("").join(seq_fragments)
            if features_sorted[0][3] == "-":
                full_seq = full_seq.reverse_complement()
            all_records.append(SeqRecord(full_seq, id=gene_id, description="CDS"))

        # Write all CDS sequences to a single FASTA
        SeqIO.write(all_records, output[0], "fasta")
