import os
from subprocess import run

configfile: "config.yaml"
groups = ["target", "base"]
GENOME= ["target_base"]

rule all:
    input:
        expand(config["genome_group"], group=groups),
        expand(config["proteome_group"], group=groups),
        expand(config["annotation_group"], group=groups),
        expand(config["mrna_group"], group=groups),

        expand(config["filtered_group"], group=groups),
        expand(config["bed_group"], group=groups),

        expand(config["annotation_db_group"], group=groups),

        expand(config["protein_db_group"], group=groups),
        config["gmap_marker"],
        config["target_cds"],
        config["base_cds"],
        config["flash_fry_index"],

        expand(config["group_pairs"], genome=GENOME),
        config["combo_gff"],
        config["genome_fna_fai"],

rule download_and_extract_data:
    output:
        genome=config["genome_group"],
        proteins=config["proteome_group"],
        annotation=config["annotation_group"],
        mrna=config["mrna_group"]
    params:
        taxon_id=lambda wc: config["taxon_ids"][wc.group]
    threads: 1
    run:
        import os, shutil, subprocess, time
        from Bio import SeqIO

        g = wildcards.group
        taxon = params.taxon_id

        gdir = f"data/genome/{g}"
        pdir = f"data/proteome/{g}"
        os.makedirs(gdir, exist_ok=True)
        os.makedirs(pdir, exist_ok=True)

        if all(os.path.exists(f) for f in output):
            print(f"[{g}] Already downloaded, skipping")
            return

        zip_tmp = f"{gdir}/bundle.zip.tmp"
        zip_final = f"{gdir}/bundle.zip"

        for f in [zip_tmp, zip_final]:
            if os.path.exists(f):
                os.remove(f)

        print(f"[{g}] Downloading taxid {taxon}")
        success = False
        for attempt in range(5):
            try:
                subprocess.run(
                    [
                        "datasets", "download", "genome", "taxon", str(taxon),
                        "--reference",
                        "--include", "genome,protein,gff3",
                        "--no-progressbar",
                        "--filename", zip_tmp
                    ],
                    check=True
                )
                os.rename(zip_tmp, zip_final)
                subprocess.run(["unzip", "-tq", zip_final], check=True)
                success = True
                break
            except subprocess.CalledProcessError:
                print(f"[{g}] Download failed (attempt {attempt+1}), retrying...")
                for f in [zip_tmp, zip_final]:
                    if os.path.exists(f):
                        os.remove(f)
                time.sleep(10)

        if not success:
            raise RuntimeError(f"[{g}] Failed to download taxid {taxon}")

        # Extract
        subprocess.run(["unzip", "-oq", zip_final, "-d", gdir], check=True)

        # Locate files deterministically
        def find_one(pattern):
            for root, _, files in os.walk(gdir):
                for f in files:
                    if f.endswith(pattern):
                        return os.path.join(root, f)
            raise FileNotFoundError(pattern)

        shutil.move(find_one("_genomic.fna"), output.genome)
        shutil.move(find_one("protein.faa"), output.proteins)
        shutil.move(find_one(".gff"), output.annotation)

        # Cleanup
        shutil.rmtree(f"{gdir}/ncbi_dataset", ignore_errors=True)
        os.remove(zip_final)

        # --- mRNA extraction directly from CDS ---
        genome_dict = SeqIO.to_dict(SeqIO.parse(output.genome, "fasta"))

        protein_cds_seqs = {}

        with open(output.annotation) as gff_file:
            for line in gff_file:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                chrom, source, ftype, start, end, score, strand, phase, attr = parts
                if ftype != "CDS":
                    continue
                attrs = {kv.split("=")[0]: kv.split("=")[1] for kv in attr.split(";") if "=" in kv}
                if "protein_id" not in attrs:
                    continue
                pid = attrs["protein_id"]
                seq = genome_dict[chrom].seq[int(start)-1:int(end)]
                if strand == "-":
                    seq = seq.reverse_complement()
                protein_cds_seqs.setdefault(pid, []).append(str(seq))

        # Concatenate CDS segments
        for pid in protein_cds_seqs:
            protein_cds_seqs[pid] = "".join(protein_cds_seqs[pid]).upper()

        # Write FASTA with headers matching protein IDs
        with open(output.mrna, "w") as out_f:
            for pid, seq in protein_cds_seqs.items():
                out_f.write(f">{pid}\n")
                for i in range(0, len(seq), 80):
                    out_f.write(seq[i:i+80] + "\n")

        if not os.path.exists(output.mrna) or os.path.getsize(output.mrna) == 0:
            raise RuntimeError(f"[{g}] mRNA FASTA generation failed")

rule filter_gff:
    input: config["annotation_group"]
    output: config["filtered_group"]
    shell:
        "awk '$3==\"gene\" || $3==\"mRNA\" || $3==\"CDS\"' {input} > {output}"
rule gff_to_bed:
    input: config["filtered_group"]
    output: config["bed_group"]
    shell: r'''
        awk '$3=="gene" {{
            split($9,a,";");
            for(i in a) if(a[i] ~ /^ID=/) {{ sub(/^ID=/,"",a[i]); id=a[i]; }}
            print $1"\t"$4-1"\t"$5"\t"id"\t.\t"$7
        }}' {input} > {output}
    '''
rule build_gff_db:
    input: config["filtered_group"]
    output: config["annotation_db_group"]
    threads: 1
    run:
        import gffutils
        gff_file = input[0]
        db_file = output[0]

        if os.path.exists(db_file):
            print(f"[{wildcards.group}] gffutils DB exists, skipping")
        else:
            print(f"[{wildcards.group}] Building gffutils DB")
            gffutils.create_db(
                gff_file,
                dbfn=db_file,
                force=True,
                keep_order=True,
                merge_strategy="merge",
                sort_attribute_values=True
            )
rule gmap_build_db:
    input: 
        genome = config["target_genome"]
    output: 
        done = config["gmap_marker"]
    threads: 4
    shell:
        """
        mkdir -p data/gmap_db/target
        # Use the conda-installed gmap_build
        gmap_build -D data/gmap_db/target -d db {input.genome}
        touch {output.done}
        """
rule make_proteome_db:
    input: config["proteome_group"]
    output: config["protein_db_group"]
    threads: 4
    shell:
        """
        makeblastdb -in {input} -dbtype prot -out data/proteome/{wildcards.group}/proteins
        touch {output}
        """

rule extract_all_cds:
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

rule flashfry_index:
    input:
        ref = config["target_genome"]
    output:
        idx = config["flash_fry_index"],
        header = config["header"]
    shell:
        """
        mkdir -p data/flashfry_tmp
        mkdir -p data/flashfry_index
        
        java -Xmx4g -jar preset/FlashFry.jar index \
            --reference {input.ref} \
            --database data/flashfry_index/genome \
            --tmpLocation data/flashfry_tmp \
            --enzyme spcas9ngg
        """

rule clean_fasta:
    input:
        fasta=lambda wc: f"data/proteome/{wc.group}/proteins.faa",
        gff=lambda wc: f"data/genome/{wc.group}/annotation.filtered.gff3"
    output:
        clean="data/mcscanx/{group}.clean.faa"
    run:
        import os
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.clean), exist_ok=True)

        # Step 1: Build protein → gene mapping from GFF
        prot_to_gene = {}
        with open(input.gff) as gff_file:
            for line in gff_file:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attrs = dict(x.split("=", 1) for x in parts[8].split(";") if "=" in x)
                gene_id = attrs.get("ID") or attrs.get("Name") or attrs.get("gene")
                parent  = attrs.get("Parent")
                if gene_id and parent:
                    prot_to_gene[parent] = gene_id

        print(f"[{wildcards.group}] Found {len(prot_to_gene)} protein→gene mappings")

        # Step 2: Rewrite FASTA headers
        written = 0
        with open(output.clean, "w") as out_f:
            for rec in SeqIO.parse(input.fasta, "fasta"):
                new_id = prot_to_gene.get(rec.id, rec.id)
                rec.id = new_id
                rec.description = ""
                SeqIO.write(rec, out_f, "fasta")
                written += 1

        if written == 0:
            raise ValueError(
                f"[{wildcards.group}] clean_fasta produced ZERO sequences — "
                "check FASTA/GFF ID compatibility"
            )

        print(f"[{wildcards.group}] Wrote {written} sequences to {output.clean}")
rule parse_gff:
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
            raise ValueError(f"[{wildcards.group}] parse_gff produced ZERO entries — no CDS protein_id/Name found")

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

        print(f"[{wildcards.group}] Wrote {len(gene_coords)} collapsed CDS entries to {output.parsed}")

rule merge_gff:
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
rule mcscanx:
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
rule diamond_make_db:
    input:
        fasta=config["clean_faa"]
    output:
        db=config["blast_dmnd"]
    threads: 4
    shell:
        """
        mkdir -p data/blast
        if [ ! -f {output.db} ]; then
            diamond makedb --in {input.fasta} -d {output.db}
        fi
        """

rule diamond_blastp:
    input:
        query=config["query_clean"],
        db=config["blast_db"]
    output:
        out=config["output_diamond"]
    threads: 8
    params:
        evalue=1e-3,
        max_target_seqs=500,
        sensitive="--more-sensitive"
    shell:
        """
        mkdir -p data/blast
        diamond blastp \
            --query {input.query} \
            --db {input.db} \
            --out {output.out} \
            --evalue {params.evalue} \
            --max-target-seqs {params.max_target_seqs} \
            --outfmt 6 \
            {params.sensitive} \
            --threads {threads}
        """
rule combine_blast:
    input:
        blast1=config["target_vs_base"],
        blast2=config["base_vs_target"]
    output:
        combined=config["combo_blast"]
    shell:
        """
        mkdir -p data/mcscanx
        cat {input.blast1} {input.blast2} > {output.combined}
        """

rule extract_all_pairs:
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
                "[all] extract_all_pairs found ZERO pairs — "
                "parser did not detect any alignment lines"
            )

        with open(output.pairs, "w") as out:
            out.write("gene1\tgene2\n")
            for g1, g2 in pairs:
                out.write(f"{g1}\t{g2}\n")

        print(f"[all] Extracted {len(pairs)} syntenic gene pairs")
rule extract_wgd_pairs:
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

rule samtools_index:
    input:
        config["target_genome"]
    output:
        config["genome_fna_fai"]
    shell:
        "samtools faidx {input}"
