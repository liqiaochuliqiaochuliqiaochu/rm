import os
from Bio import SeqIO

# Constants
flank = 1000 #change this value to change upstream/downstream length
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

        # Load GenBank file
        try:
            record = SeqIO.read(gbk_path, "genbank")
        except:
            record = next(SeqIO.parse(gbk_path, "genbank"))
            print(f"Warning: multiple genbnk records exist for {base_name}, please manually check." )
        genome_seq = record.seq
        genome_len = len(genome_seq)

        with open(result_path) as infile, open(output_path, "w") as outfile:
            for line in infile:
                full_line = line.strip()
                if not full_line:
                    continue
                try:
                    # Get position from the first word in the line
                    loc = full_line.split()[0]
                    start_str, end_str = loc.split("_")
                    start, end = int(start_str), int(end_str)

                    # Calculate sequence boundaries
                    region_start = max(0, start - flank)
                    region_end = min(genome_len, end + flank)

                    # Extract gene and flanking sequences
                    gene_seq = genome_seq[start:end]
                    region_seq = genome_seq[region_start:region_end]

                    # Detect strand direction
                    strand = None
                    for feature in record.features:
                        if feature.type == "gene":
                            if (feature.location.start == start and
                                feature.location.end == end):
                                strand = feature.location.strand
                                break

                    # Reverse complement if on negative strand
                    if strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                        region_seq = region_seq.reverse_complement()

                    # Write gene-only sequence
                    outfile.write(f">{full_line}\n{gene_seq}\n")
                    # Write sequence with flanks
                    outfile.write(f">{full_line} +flank\n{region_seq}\n")

                except Exception as e:
                    print(f"Error processing line: {full_line} in {result_path}")
                    print(e)
