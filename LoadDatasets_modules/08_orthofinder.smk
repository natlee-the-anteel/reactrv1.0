# LoadDatasets_modules/08_orthofinder.smk
# Purpose: Identify orthologous proteins using OrthoFinder
# Note: OrthoFinder is optional - can be skipped if not needed

rule orthofinder:
    """Detect orthologous proteins between target and base genomes.
    
    OrthoFinder performs whole-genome ortholog identification using
    bidirectional BLAST search and graph-based clustering.
    
    Inputs: Protein sequences from both genomes
    Output: OrthoFinder results directory with ortholog groups (Orthogroups.tsv)
    
    Parameters:
    - method: DIAMOND (fast) or BLAST (more compatible)
    - threads: Number of parallel DIAMOND/BLAST searches
    
    Note: Run this separately or include in full pipeline as needed:
        snakemake -s LoadDatasets.smk -R orthofinder extract_orthogroups -j 8
    """
    input:
        target_prot = config["target_proteome"],
        base_prot = config["base_proteome"]
    output:
        ortho_marker = "data/blast/orthofinder.done"
    params:
        method = "diamond",
        threads = 8
    shell:
        """
        mkdir -p data/blast
        
        # Run OrthoFinder on the proteome directory
        # -f: input directory (data/proteome contains target/ and base/ subdirs with *.faa)
        # -t: number of threads
        # -M: sequence search method (diamond or blast)
        # -o: output directory
        cd data/blast
        orthofinder -f ../proteome \
            -t {params.threads} \
            -M {params.method} \
            -o ./OrthoFinder \
            -op
        
        cd ../..
        touch data/blast/orthofinder.done
        """

rule extract_orthogroups:
    """Extract single-copy orthologous gene pairs from OrthoFinder results.
    
    Filters OrthoFinder output to keep only 1-to-1 orthologs
    (single-copy genes present in both genomes exactly once).
    
    Output: TSV file with columns: target_gene, base_gene, orthogroup_id
    """
    input:
        ortho_marker = "data/blast/orthofinder.done",
        orthogroups = "data/blast/OrthoFinder/Results/Orthogroups.tsv"
    output:
        ortho_pairs = "data/blast/orthofinder_1to1_pairs.tsv"
    run:
        import os
        import pandas as pd

        # Read OrthoFinder's Orthogroups.tsv
        # Columns: Orthogroup, Species1, Species2, ...
        df = pd.read_csv(input.orthogroups, sep="\t")
        
        # Rename columns to make them more generic (handles multiple species)
        cols = df.columns.tolist()
        species_cols = [c for c in cols if c != "Orthogroup"]
        
        os.makedirs(os.path.dirname(output.ortho_pairs), exist_ok=True)
        
        # Extract 1-to-1 orthologs
        one_to_one = []
        for _, row in df.iterrows():
            # Skip orthogroups that don't have exactly 2 species
            if len(species_cols) < 2:
                continue
            
            species1_genes = str(row[species_cols[0]]).split(", ") if pd.notna(row[species_cols[0]]) else []
            species2_genes = str(row[species_cols[1]]).split(", ") if pd.notna(row[species_cols[1]]) else []
            
            # Keep only if 1 gene in each species (1-to-1 ortholog)
            if len(species1_genes) == 1 and len(species2_genes) == 1 and species1_genes[0] != "" and species2_genes[0] != "":
                one_to_one.append({
                    "target_gene": species1_genes[0],
                    "base_gene": species2_genes[0],
                    "orthogroup_id": row["Orthogroup"]
                })
        
        # Write output
        with open(output.ortho_pairs, "w") as out:
            out.write("target_gene\tbase_gene\torthogroup_id\n")
            for pair in one_to_one:
                out.write(f"{pair['target_gene']}\t{pair['base_gene']}\t{pair['orthogroup_id']}\n")
        
        print(f"Extracted {len(one_to_one)} 1-to-1 ortholog pairs")
