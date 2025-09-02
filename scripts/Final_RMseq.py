import os
import glob
from Bio import SeqIO
import re
import csv

# --- Config ---
rm_folder = "FinalRM"
gbk_folder = "genbank"
output_folder = "sequence"
list_file = "cleanRM.txt"
merged_file = "RMall.fasta"
summary_file = "RM_summary.csv"
type_count_file = "RM_type_counts.csv"
CLUSTER_WINDOW = 10000  # 10 kb max distance for clustering

os.makedirs(output_folder, exist_ok=True)

def clean_sequence(seq):
    return ''.join([base for base in str(seq).upper() if base in "ATCGN"])

def parse_coordinate(coord_str):
    match = re.search(r'\d+', coord_str)
    return int(match.group()) if match else None

# --- Main ---
with open(merged_file, "w") as merged_out, \
     open(summary_file, "w", newline="") as sum_out, \
     open(type_count_file, "w", newline="") as type_out:

    sum_writer = csv.writer(sum_out)
    sum_writer.writerow(["Genome", "TypeI_complete", "TypeII", "TypeIII", "TypeIV", "Total_RM"])

    type_writer = csv.writer(type_out)
    type_writer.writerow([
        "Genome", "IR", "IM", "IS", "IRM",
        "IIR", "IIRM", "IIGR", "IIGM", "IIIM", "IIIR", "IIIRM", "IVR"
    ])

    with open(list_file) as lf:
        for line in lf:
            sample_full = line.strip()
            if not sample_full:
                continue
            print(f"Processing: {sample_full}")

            base_name = sample_full.split("\t")[0].replace("_genomic", "")

            # Result file
            result_files = glob.glob(os.path.join(rm_folder, base_name + "*.txt"))
            if not result_files:
                print(f"Result file not found for {base_name}, skipping.")
                continue
            result_path = result_files[0]

            # GenBank file
            gbk_files = glob.glob(os.path.join(gbk_folder, base_name + "*.gbk")) + \
                        glob.glob(os.path.join(gbk_folder, base_name + "*.gbff"))
            if not gbk_files:
                print(f"GenBank file not found for {sample_full}, skipping.")
                continue
            gbk_path = gbk_files[0]

            # Load records
            records = list(SeqIO.parse(gbk_path, "genbank"))
            if not records:
                print(f"No sequences found in {gbk_path}, skipping.")
                continue

            output_path = os.path.join(output_folder, base_name + ".fasta")

            typeI_genes, typeII_genes, typeIII_genes, typeIV_genes = [], [], [], []

            # Track counts for each type
            counts = {k: 0 for k in ["IR","IM","IS","IRM","IIR","IIRM","IIGR","IIGM","IIIM","IIIR","IIIRM","IVR"]}

            # Read RM prediction file
            with open(result_path) as infile:
                for line in infile:
                    parts = line.strip().split()
                    if len(parts) < 3:
                        continue
                    gene_type = parts[2]
                    loc = parts[0]
                    coords = loc.split("_")
                    start, end = map(parse_coordinate, coords[-2:])
                    if start is None or end is None:
                        continue

                    if gene_type in ["IR","IM","IS"]:
                        typeI_genes.append((gene_type,start,end,parts))
                    elif gene_type in ["IIIR","IIIM","IIIRM"]:
                        typeIII_genes.append((gene_type,start,end,parts))
                    elif gene_type in ["IIR","IIRM","IIGR","IIGM"]:
                        typeII_genes.append((gene_type,start,end,parts))
                    elif gene_type == "IVR":
                        typeIV_genes.append((gene_type,start,end,parts))

                    if gene_type in counts:
                        counts[gene_type] += 1

            # --- Type I: sequential + bidirectional scan ---
            TypeI_complete_count = 0
            complete_typeI_R = []
            used_genes = set()

            for i, g1 in enumerate(typeI_genes):
                if (g1[1], g1[2]) in used_genes or g1[0] != "IR":
                    continue

                candidate = [g1]
                used_genes.add((g1[1], g1[2]))
                found_M = None
                found_S = None

                # Scan downstream
                for g2 in typeI_genes[i+1:]:
                    if (g2[1], g2[2]) in used_genes:
                        continue
                    if g2[1] - g1[2] > CLUSTER_WINDOW:
                        break
                    if g2[0] == "IM" and not found_M:
                        found_M = g2
                        candidate.append(g2)
                        used_genes.add((g2[1], g2[2]))
                    elif g2[0] == "IS" and not found_S:
                        found_S = g2
                        candidate.append(g2)
                        used_genes.add((g2[1], g2[2]))
                    elif g2[0] == "IR":
                        break

                # Scan upstream
                for g2 in reversed(typeI_genes[:i]):
                    if (g2[1], g2[2]) in used_genes:
                        continue
                    if g1[1] - g2[2] > CLUSTER_WINDOW:
                        break
                    if g2[0] == "IM" and not found_M:
                        found_M = g2
                        candidate.append(g2)
                        used_genes.add((g2[1], g2[2]))
                    elif g2[0] == "IS" and not found_S:
                        found_S = g2
                        candidate.append(g2)
                        used_genes.add((g2[1], g2[2]))

                if found_M and found_S:
                    TypeI_complete_count += 1
                    complete_typeI_R.append(g1)

            # --- Type III: IIIR + IIIM bidirectional OR standalone IIIRM ---
            TypeIII_complete_R = []
            TypeIII_complete_count = 0
            used_genes_III = set()
            typeIII_genes.sort(key=lambda x: x[1])

            for i, g1 in enumerate(typeIII_genes):
                if g1[0] != "IIIR" or (g1[1], g1[2]) in used_genes_III:
                    continue
                found_M = None

                # Scan downstream
                for g2 in typeIII_genes[i+1:]:
                    if (g2[1], g2[2]) in used_genes_III:
                        continue
                    if g2[1] - g1[2] > CLUSTER_WINDOW:
                        break
                    if g2[0] == "IIIM":
                        found_M = g2
                        used_genes_III.add((g2[1], g2[2]))
                        break

                # Scan upstream
                if not found_M:
                    for g2 in reversed(typeIII_genes[:i]):
                        if (g2[1], g2[2]) in used_genes_III:
                            continue
                        if g1[1] - g2[2] > CLUSTER_WINDOW:
                            break
                        if g2[0] == "IIIM":
                            found_M = g2
                            used_genes_III.add((g2[1], g2[2]))
                            break

                if found_M:
                    TypeIII_complete_count += 1
                    TypeIII_complete_R.append(g1)
                    used_genes_III.add((g1[1], g1[2]))

            # Standalone IIIRM
            for g in typeIII_genes:
                if g[0] == "IIIRM" and (g[1], g[2]) not in used_genes_III:
                    TypeIII_complete_count += 1
                    TypeIII_complete_R.append(g)
                    used_genes_III.add((g[1], g[2]))

            # --- Type II: independent (IIR, IIRM, IIGR, IIGM) ---
            TypeII_count = counts["IIR"] + counts["IIRM"] + counts["IIGR"] + counts.get("IIGM",0)
            TypeII_R = [g for g in typeII_genes if g[0] in ["IIR","IIRM","IIGR","IIGM"]]

            # --- Type IV ---
            TypeIV_count = counts["IVR"]
            TypeIV_R = typeIV_genes[:]

            # Collect genes to export
            export_genes = complete_typeI_R + TypeII_R + TypeIII_complete_R + TypeIV_R

            # --- Extract sequences ---
            with open(output_path,"w") as outfile:
                for gene in export_genes:
                    gene_type, start, end, parts = gene
                    gene_seq = None
                    for record in records:
                        if start <= len(record.seq) and end <= len(record.seq):
                            gene_seq = record.seq[start-1:end]
                            strand = None
                            for feature in record.features:
                                if feature.type=="gene":
                                    f_start=int(feature.location.start)+1
                                    f_end=int(feature.location.end)
                                    if f_start==start and f_end==end:
                                        strand = feature.location.strand
                                        break
                            if strand == -1:
                                gene_seq = gene_seq.reverse_complement()
                            break
                    if gene_seq is None:
                        print(f"Gene {start}-{end} not found in {sample_full}")
                        continue
                    clean_seq = clean_sequence(gene_seq)
                    header_str = re.sub(r'[^A-Za-z0-9_]','_',f"{base_name}_{parts[1]}_{gene_type}")
                    outfile.write(f">{header_str}\n{clean_seq}\n")
                    merged_out.write(f">{header_str}\n{clean_seq}\n")

            # --- Summaries ---
            total_rm = TypeI_complete_count + TypeII_count + TypeIII_complete_count + TypeIV_count
            sum_writer.writerow([
                base_name,
                TypeI_complete_count,
                TypeII_count,
                TypeIII_complete_count,
                TypeIV_count,
                total_rm
            ])
            type_writer.writerow([
                base_name,
                counts["IR"], counts["IM"], counts["IS"], counts["IRM"],
                counts["IIR"], counts["IIRM"], counts["IIGR"], counts.get("IIGM",0),
                counts.get("IIIM",0), counts.get("IIIR",0), counts.get("IIIRM",0), counts["IVR"]
            ])

print("Processing complete. RM sequences extracted, summary and type counts exported.")
