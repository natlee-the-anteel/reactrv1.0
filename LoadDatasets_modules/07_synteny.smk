# LoadDatasets_modules/07_synteny.smk
# Purpose: Identify syntenic blocks and classify gene duplicates

rule mcscanx:
    """Run MCScanX to identify syntenic blocks.
    
    Parameters:
    - s=2: Number of sequences (always 2 for target vs base)
    - g=500: Gap penalty (larger = allow more distant genes in block)
    - m=10: Minimum genes per block
    - e=0.01: E-value threshold
    
    Output: collinearity file describing block structure and gene pairs.
    """
    input:
        gff=config["group_filtered_gff"],
        blast=config["group_blast"]
    output:
        collinearity=config["group_collin"]
    params:
        s=2,
        g=500,
        m=10,
        e=0.01
    threads: 8
    shell:
        """
        preset/MCScanX data/mcscanx/{wildcards.genome} \
            -s {params.s} -g {params.g} -m {params.m} -e {params.e}
        """

rule duplicate_gene_classifier:
    """Classify genes as duplicated or singleton based on MCScanX output."""
    input:
        collinearity=config["group_collin"],
        gff=config["group_filtered_gff"]
    output:
        gene_type=config["gene_type"]
    threads: 1
    shell:
        """
        preset/duplicate_gene_classifier data/mcscanx/{wildcards.genome}
        """

rule extract_all_pairs:
    """Parse MCScanX collinearity file to extract all syntenic gene pairs.
    
    Reads collinearity file and extracts all (gene1, gene2) pairs
    from alignment blocks. Output: TSV with gene1 and gene2 columns.
    """
    input:
        collinearity=config["combo_collin"]
    output:
        pairs=config["all_pair"]
    run:
        import os

        os.makedirs(os.path.dirname(output.pairs), exist_ok=True)
        pairs = []

        with open(input.collinearity) as f:
            for line in f:
                line = line.strip()

                # Skip headers and alignment titles
                if (
                    not line
                    or line.startswith("#")
                    or line.startswith("## Alignment")
                ):
                    continue

                # Valid gene-pair lines always contain ":"
                # Example:
                # 0-  2: XXG40820.1 NP_001323410.1 0.0005
                if ":" not in line:
                    continue

                parts = line.split()

                # Defensive check
                if len(parts) < 4:
                    continue

                gene1 = parts[2]
                gene2 = parts[3]

                pairs.append((gene1, gene2))

        if not pairs:
            raise ValueError(
                "extract_all_pairs found ZERO pairs — "
                "parser did not detect any alignment lines"
            )

        with open(output.pairs, "w") as out:
            out.write("gene1\tgene2\n")
            for g1, g2 in pairs:
                out.write(f"{g1}\t{g2}\n")

        print(f"Extracted {len(pairs)} syntenic gene pairs")

rule extract_wgd_pairs:
    """Extract whole-genome duplication (WGD) pairs from collinearity.
    
    Identifies genes duplicated by WGD events (large-scale duplications)
    and classifies them by duplication block.
    """
    input:
        collinearity=config["group_collin"]
    output:
        pairs=config["group_pairs"]
    shell:
        r"""
        awk '
        /^## Alignment/ {{
            block = $3
            sub(/:$$/, "", block)
        }}
        !/^#/ && NF >= 4 {{
            print $3 "\t" $4 "\t" block
        }}
        ' {input.collinearity} > {output.pairs}
        """
