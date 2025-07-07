import os
from Bio import SeqIO

# Constants
flank = 2000  # change this value to change upstream/downstream length
rm_folder = "FinalRM"
gbk_folder = "gbk"
output_folder = "sequence"
list_file = "list.txt"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Read names from list.txt
with open(list_file) as lf:
    for line in lf:
        if not line.strip():
            continue
        base_name = line.strip().split(".")[0]
        print(base_name)
        result_path = os.path.join(rm_folder, base_name + ".txt")
        gbk_path = os.path.join(gbk_folder, base_name + ".gbk")
        output_path = os.path.join(output_folder, base_name + ".fasta")

        if not os.path.exists(result_path):
            print(f"Result file not found for {base_name}, skipping.")
            continue
        if not os.path.exists(gbk_path):
            print(f"GenBank file not found for {base_name}, skipping.")
            continue

        try:
            records = list(SeqIO.parse(gbk_path, "genbank"))
        except Exception as e:
            print(f"Error reading {gbk_path}: {e}")
            continue

        with open(result_path) as infile, open(output_path, "w") as outfile:
            for line in infile:
                full_line = line.strip()
                if not full_line:
                    continue
                try:
                    loc = full_line.split()[0]
                    start_str, end_str = loc.split("_")
                    start, end = int(start_str), int(end_str)
                    gene_found = False

                    for record in records:
                        for feature in record.features:
                            if feature.type == "CDS":
                                feat_start = int(feature.location.start)
                                feat_end = int(feature.location.end)

                                # Match by position
                                if feat_start == start - 1 and feat_end == end:
                                    gene_seq = feature.extract(record.seq)

                                    # Determine strand
                                    strand = feature.location.strand
                                    if strand == -1:
                                        gene_seq = gene_seq.reverse_complement()

                                    # Flanking region (based on genomic positions)
                                    genome_len = len(record.seq)
                                    region_start = max(0, start - flank)
                                    region_end = min(genome_len, end + flank)
                                    region_seq = record.seq[region_start:region_end]
                                    if strand == -1:
                                        region_seq = region_seq.reverse_complement()

                                    # Write outputs
                                    outfile.write(f">{full_line}\n{gene_seq}\n")
                                    outfile.write(f">{full_line} +flank\n{region_seq}\n")

                                    gene_found = True
                                    break
                        if gene_found:
                            break

                    if not gene_found:
                        print(f"CDS {start}_{end} not found in any contig for {base_name}")

                except Exception as e:
                    print(f"Error processing line: {full_line} in {result_path}")
                    print(e)
