from Bio import Seq


def read_content(filename):
    """
    Reads a tab-delimited file and returns its contents as a list of lists.
    
    :param filename: Path to the file to be read.
    :return: A list of lists, where each sublist represents a line in the file.
    """
    with open(filename, 'r') as file:
        return [line.strip().split('\t') for line in file]


# Read input files
list_data = read_content("list.txt")
ref = read_content("/home5/qli37/24RM/final/dna.txt")

# Open output file for writing results
with open("FinalRM.txt", 'w') as f1:
    f1.write("name\tTypeI\tTypeII\tTypeIII\tTypeIV\tTotal\n")

    for entry in list_data:
        try:
            name = entry[0].split(".gbk")[0]
            arr = []
            seq_name = set()

            protein_result = read_content(f"BLASTP/{name}_clean.blast")
            ip = read_content(f"aaIn/{name}_clean.fasta")
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue

        for row in protein_result:
            try:
                if ((row[0] not in seq_name) and float(row[2]) > 75):  
                    t, st = "Unknown", ""

                    # Find RM system type from the reference file
                    for i in range(0, len(ref) - 1, 2):
                        info = ref[i][0].split(">REBASE:")[1]
                        ref_name = info.split("_")[0]

                        if row[1].split("_")[1] == ref_name:
                            try:
                                t = info.split("Type_")[1].split("_")[0]
                            except:
                                t = 'orphan' if "orphan" in info else "Unknown"

                            # Determine subtype
                            if "/" in ref_name:
                                st = "RM"
                            elif "." in ref_name:
                                st = ref_name.split(".")[0]
                                st = 'M' if 'M' in st else 'S' if 'S' in st else 'R' if 'R' in st else 'C'
                                if st == 'C':
                                    t = "III"
                            elif "restriction" in info:
                                st = "R"
                            elif "methyltransferase" in info:
                                st = "M"
                            elif "control" in info:
                                st, t = "C", "III"
                            break

                    arr.append([row[0], row[1], t, st, "_".join(row[2:6])])
                    seq_name.add(row[0])

            except Exception as e:
                print(f"Error processing row {row}: {e}")
                continue

        # Write individual results
        with open(f"FinalRM/{name}.txt", 'w') as f:
            type_counts = {
                "IM": 0, "IS": 0, "IR": 0, "IRM": 0, "IIM": 0, "IIR": 0, "IIRM": 0,
                "IIS": 0, "IIIM": 0, "IIIS": 0, "IIIR": 0, "IIIRM": 0, "IVR": 0,
                "Other": 0, "orphan1": 0, "orphan2": 0, "orphan3": 0, "orphan4": 0
            }

            for record in arr:
                f.write("\t".join(record) + "\n")
                ty = record[2] + record[3]

                if ty in type_counts:
                    type_counts[ty] += 1
                else:
                    type_counts["Other"] += 1

            # Calculate type classifications
            typeI = min(type_counts["IR"], type_counts["IM"], type_counts["IS"]) if type_counts["IR"] else type_counts["orphan1"]
            typeII = type_counts["IIR"] if type_counts["IIR"] else type_counts["orphan2"]
            typeIII = min(type_counts["IIIR"], type_counts["IIIM"]) if type_counts["IIIR"] else type_counts["orphan3"]
            typeIV = type_counts["IVR"] if type_counts["IVR"] else type_counts["orphan4"]

            total = typeI + type_counts["IRM"] + typeII + type_counts["IIRM"] + typeIII + type_counts["IIIRM"] + typeIV

            # Write summary to FinalRM.txt
            f1.write(f"{name}\t{typeI}\t{typeII}\t{typeIII}\t{typeIV}\t{total}\n")
