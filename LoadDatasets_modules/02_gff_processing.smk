# LoadDatasets_modules/02_gff_processing.smk
# Purpose: GFF3 filtering, format conversion, and preparation for MCScanX

rule filter_gff:
    """Filter GFF3 to keep only gene, mRNA, and CDS features."""
    input: config["annotation_group"]
    output: config["filtered_group"]
    shell:
        "awk '$3==\"gene\" || $3==\"mRNA\" || $3==\"CDS\"' {input} > {output}"

rule gff_to_bed:
    """Convert filtered GFF3 to BED format (gene features only)."""
    input: config["filtered_group"]
    output: config["bed_group"]
    shell: r'''
        awk '$3=="gene" {{
            split($9,a,";");
            for(i in a) if(a[i] ~ /^ID=/) {{ sub(/^ID=/,"",a[i]); id=a[i]; }}
            print $1"\t"$4-1"\t"$5"\t"id"\t.\t"$7
        }}' {input} > {output}
    '''

rule parse_gff:
    """Parse GFF3 to extract CDS coordinates for MCScanX format.
    
    Collapses multiple CDS per gene into single min/max coordinates.
    Output: chromosome, gene_id, start, end (tab-separated).
    """
    input:
        gff=config["filtered_group"]
    output:
        parsed=config["scan_combo"]
    run:
        import os

        os.makedirs(os.path.dirname(output.parsed), exist_ok=True)
        rows = []

        with open(input.gff) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                chrom = parts[0]
                start, end = int(parts[3]), int(parts[4])
                attrs = dict(x.split("=", 1) for x in parts[8].split(";") if "=" in x)
                gene_id = attrs.get("protein_id") or attrs.get("Name")
                if gene_id:
                    rows.append((chrom, gene_id, start, end))

        if not rows:
            raise ValueError(f"parse_gff produced ZERO entries — no CDS protein_id/Name found")

        # Collapse multiple CDS per gene
        gene_coords = {}
        for chrom, gene_id, start, end in rows:
            if gene_id in gene_coords:
                gene_coords[gene_id]["start"] = min(gene_coords[gene_id]["start"], start)
                gene_coords[gene_id]["end"] = max(gene_coords[gene_id]["end"], end)
            else:
                gene_coords[gene_id] = {"chrom": chrom, "start": start, "end": end}

        with open(output.parsed, "w") as out:
            for gene_id, vals in gene_coords.items():
                out.write(f"{vals['chrom']}\t{gene_id}\t{vals['start']}\t{vals['end']}\n")

        print(f"Wrote {len(gene_coords)} collapsed CDS entries to {output.parsed}")

rule merge_gff:
    """Merge target and base GFF files for combined MCScanX analysis."""
    input:
        target=config["filtered_target_gff"],
        base=config["filtered_base_gff"]
    output:
        merged=config["combo_gff"]
    shell:
        """
        mkdir -p data/mcscanx
        cat {input.target} {input.base} > {output.merged}
        """
