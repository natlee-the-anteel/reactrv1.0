# LoadDatasets_modules/01_download.smk
# Purpose: Download and extract genomes, proteomes, and annotations from NCBI

rule download_and_extract_data:
    """Download genome, proteome, and annotation from NCBI datasets API.
    
    Includes robust retry logic (5 attempts) and post-processing:
    - GFF3 sorting by chromosome/start/type/length
    - Collapsing multi-part (trans-spliced) gene entries
    - mRNA extraction from CDS annotations
    """
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

        # Cleanup temp files
        if os.path.exists(zip_tmp): os.remove(zip_tmp)
        if os.path.exists(zip_final): os.remove(zip_final)

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
                if os.path.exists(zip_tmp): os.remove(zip_tmp)
                time.sleep(10)

        if not success:
            raise RuntimeError(f"[{g}] Failed to download taxid {taxon}")

        # Extract
        subprocess.run(["unzip", "-oq", zip_final, "-d", gdir], check=True)

        # Helper to find files
        def find_one(pattern):
            for root, _, files in os.walk(gdir):
                for f in files:
                    if f.endswith(pattern):
                        return os.path.join(root, f)
            raise FileNotFoundError(pattern)

        shutil.move(find_one("_genomic.fna"), output.genome)
        shutil.move(find_one("protein.faa"), output.proteins)
        
        # --- ROBUST PYTHON GFF3 SORTING ---
        # We use Python here to avoid OS-specific differences in the 'sort' command
        raw_gff = find_one(".gff")
        print(f"[{g}] Sorting GFF3 headers and features...")
        
        with open(raw_gff, 'r') as fin:
            lines = fin.readlines()

        headers = [line for line in lines if line.startswith("#")]
        body = [line for line in lines if not line.startswith("#")]

        def gff_sort_key(line):
            parts = line.split('\t')
            if len(parts) < 9: return ('', 0, 0, 0)  # Safety
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            ftype = parts[2]

            # Define explicit type priority so that `gene` always comes before
            # `mRNA`, and `mRNA` before `exon`/`CDS`. This prevents child
            # features being listed before their parents which can break
            # downstream parsers (e.g. BCBio's remap_id) for split/trans-
            # spliced genes.
            if ftype == 'gene':
                type_prio = 0
            elif ftype == 'mRNA':
                type_prio = 1
            else:
                type_prio = 2

            length = end - start

            # Sort by: chrom, start, type priority (gene->mRNA->others),
            # then by descending length to keep longer parent features before
            # shorter ones when type is equal.
            return (chrom, start, type_prio, -length)

        body.sort(key=gff_sort_key)

        # Collapse multi-part `gene` entries (e.g., trans-spliced genes)
        # into a single gene spanning min(start)..max(end). This ensures
        # parent `gene` ranges encompass their child `mRNA`/`CDS` and
        # avoids remap_id failures in downstream parsers (BCBio).
        genes_by_id = {}
        others = []
        for line in body:
            parts = line.split('\t')
            if len(parts) < 9:
                others.append(line)
                continue
            ftype = parts[2]
            attrs = parts[8]
            # extract ID if present
            gid = None
            for token in attrs.split(';'):
                if token.startswith('ID='):
                    gid = token.split('=', 1)[1]
                    break
            if ftype == 'gene' and gid:
                start = int(parts[3])
                end = int(parts[4])
                genes_by_id.setdefault(gid, []).append((line, parts, start, end, attrs))
            else:
                others.append(line)

        new_genes = []
        for gid, entries in genes_by_id.items():
            if len(entries) == 1:
                new_genes.append(entries[0][0])
            else:
                chrom = entries[0][1][0]
                src = entries[0][1][1]
                strand = entries[0][1][6]
                score = '.'
                minstart = min(e[2] for e in entries)
                maxend = max(e[3] for e in entries)
                attrs = entries[0][4]
                # remove part= attributes which describe fragment numbering
                attrs = ';'.join([a for a in attrs.split(';') if not a.startswith('part=')])
                new_line = '\t'.join([chrom, src, 'gene', str(minstart), str(maxend), score, str(strand), '.', attrs])
                new_genes.append(new_line)

        new_body = new_genes + others

        with open(output.annotation, 'w') as fout:
            fout.writelines(headers)
            fout.writelines([l + "\n" if not l.endswith("\n") else l for l in new_body])
        # ---------------------------

        # Cleanup download artifacts
        shutil.rmtree(f"{gdir}/ncbi_dataset", ignore_errors=True)
        os.remove(zip_final)

        # --- mRNA extraction ---
        genome_dict = SeqIO.to_dict(SeqIO.parse(output.genome, "fasta"))
        protein_cds_seqs = {}

        with open(output.annotation) as gff_file:
            for line in gff_file:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9: continue
                chrom, source, ftype, start, end, score, strand, phase, attr = parts
                if ftype != "CDS": continue
                attrs = {kv.split("=")[0]: kv.split("=")[1] for kv in attr.split(";") if "=" in kv}
                if "protein_id" not in attrs: continue
                pid = attrs["protein_id"]
                seq = genome_dict[chrom].seq[int(start)-1:int(end)]
                if strand == "-":
                    seq = seq.reverse_complement()
                protein_cds_seqs.setdefault(pid, []).append(str(seq))

        for pid in protein_cds_seqs:
            protein_cds_seqs[pid] = "".join(protein_cds_seqs[pid]).upper()

        with open(output.mrna, "w") as out_f:
            for pid, seq in protein_cds_seqs.items():
                out_f.write(f">{pid}\n")
                for i in range(0, len(seq), 80):
                    out_f.write(seq[i:i+80] + "\n")

        if not os.path.exists(output.mrna) or os.path.getsize(output.mrna) == 0:
            raise RuntimeError(f"[{g}] mRNA FASTA generation failed")
