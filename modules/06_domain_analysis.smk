# ============================================================================
# MODULE: 06_domain_analysis.smk
# Purpose: Motif analysis (MEME) and protein property calculation
# ============================================================================

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
                aa_percents = analyser.get_amino_acids_percent()
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
