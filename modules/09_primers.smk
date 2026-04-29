# ============================================================================
# MODULE: 09_primers.smk
# Purpose: Primer design and primer search analysis
# ============================================================================

rule design_primers_per_domain:
    input:
        mrna_fasta=lambda wildcards: config["mrna"].format(query=wildcards.query, domain=wildcards.domain)
    output:
        primer_file=config["primers"]
    wildcard_constraints:
        ptype="|".join(["cloning", "validation", "qpcr"])
    params:
        product_ranges={
            "qpcr":"80-150",
            "validation":"300-800",
            "cloning":"700-1500"
        }
    run:
        import os
        import shutil
        from Bio import SeqIO
        from subprocess import run, PIPE

        os.makedirs(os.path.dirname(output.primer_file), exist_ok=True)
        product_range = params.product_ranges[wildcards.ptype]

        seqs = list(SeqIO.parse(input.mrna_fasta, "fasta"))
        if not seqs:
            print(f"Warning: no mRNA sequences for domain {wildcards.domain} → primers will not be designed")
            open(output.primer_file, "w").close()
            return

        # Clear file
        open(output.primer_file, 'w').close()

        for record in seqs:
            seq_id = record.id
            sequence = str(record.seq).replace("U", "T")

            primer_input = "\n".join([
                f"SEQUENCE_ID={seq_id}",
                f"SEQUENCE_TEMPLATE={sequence}",
                "PRIMER_TASK=generic",
                "PRIMER_OPT_SIZE=20",
                "PRIMER_MIN_SIZE=18",
                "PRIMER_MAX_SIZE=25",
                "PRIMER_MIN_TM=57.0",
                "PRIMER_OPT_TM=60.0",
                "PRIMER_MAX_TM=63.0",
                f"PRIMER_PRODUCT_SIZE_RANGE={product_range}",
                "P3_FILE_FLAG=0",
                "PRIMER_EXPLAIN_FLAG=1",
                "="
            ]) + "\n"

            try:
                result = run("primer3_core", input=primer_input.encode(), stdout=PIPE, stderr=PIPE, check=False)
                with open(output.primer_file, "a") as out_f:
                    out_f.write(f"== {seq_id} ==\n")
                    out_f.write(result.stdout.decode() + "\n")
            except FileNotFoundError:
                print(f"Warning: primer3_core not found. Skipping primer design for {seq_id}")
                with open(output.primer_file, "a") as out_f:
                    out_f.write(f"== {seq_id} ==\n")
                    out_f.write("primer3_core not installed - primer design skipped\n")

rule primers_to_primersearch:
    input:
        primer3_txt=config["primers"]
    output:
        primersearch=config["primersearch_txt"]
    params:
        valid_types=["cloning", "validation", "qpcr"]
    run:
        import re

        if wildcards.ptype not in params.valid_types:
            raise ValueError(f"Invalid primer type: {wildcards.ptype}")

        out_lines = []

        with open(input.primer3_txt) as f:
            content = f.read()

        # Split by primer blocks labeled by "== geneid =="
        blocks = re.split(r"==\s*(\S+)\s*==", content)

        for i in range(1, len(blocks), 2):
            seq_id = blocks[i].strip()
            block = blocks[i + 1]

            # Find all PRIMER_LEFT_n and PRIMER_RIGHT_n entries
            lefts = re.findall(r"PRIMER_LEFT_(\d+)_SEQUENCE=([A-Za-z]+)", block)
            rights = re.findall(r"PRIMER_RIGHT_(\d+)_SEQUENCE=([A-Za-z]+)", block)

            # Pair by shared index
            pairs = {}
            for idx, seq in lefts:
                pairs[idx] = {"left": seq}
            for idx, seq in rights:
                if idx in pairs:
                    pairs[idx]["right"] = seq

            # Write output in EMBOSS primersearch format
            for idx in sorted(pairs.keys(), key=int):
                pair = pairs[idx]
                if "left" in pair and "right" in pair:
                    line = f"{seq_id}_pair{idx} {pair['left']} {pair['right']}"
                    out_lines.append(line)

        with open(output.primersearch, "w") as out_f:
            out_f.write("\n".join(out_lines) + "\n")

rule primersearch:
    input:
        genome = config["target_genome"],
        primers = config["primersearch_txt"]
    output:
        report = config["primersearch_sim"]
    params:
        mismatch = 0,       # max mismatches allowed
    shell:
        """
        primersearch -seqall {input.genome} \
                    -infile {input.primers} \
                    -mismatch {params.mismatch} \
                    -outfile {output.report}
        """
