# ============================================================================
# MODULE: 11_visualization.smk
# Purpose: Chromosome maps and gene structure visualization
# ============================================================================

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
