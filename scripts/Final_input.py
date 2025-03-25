from Bio import SeqIO

# Input and output files
def read_content(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    i = 0
    data = []
    for line in lines:
        line = line.strip().split(' ')  # strip \n, split by space
        data.append(line)  # record data
        i = i + 1
    return data

list = read_content("list.txt")
for i in range(0, len(list)):
    name = list[i][0].split(".gbk")[0]

    input_genbank = "gbk/"+name+".gbk"
    output_fasta = "In/"+name+".fasta"

    # Extract features and rename by start and end positions
    with open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(input_genbank, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":  # Change to "gene" if needed
                    # Get the start and end positions (1-based indexing)
                    start = feature.location.start + 1  # Convert to 1-based
                    end = feature.location.end

                    # Create a new name using start and end positions
                    feature_name = f"{start}_{end}"

                    # Extract the sequence
                    sequence = feature.extract(record.seq)

                    # Write to FASTA file
                    fasta_out.write(f">{feature_name}\n{sequence}\n")

    print(f"FASTA file with positions saved as {output_fasta}")

    from Bio import SeqIO

    # Input and output files

    output_fasta = "In/"+name+"_protein.fasta"

    # Extract precomputed protein sequences from the "translation" qualifier
    with open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(input_genbank, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":  # Focus on CDS features
                    # Get start and end positions for naming
                    start = feature.location.start + 1  # Convert to 1-based
                    end = feature.location.end

                    # Check for "translation" qualifier
                    if "translation" in feature.qualifiers:
                        protein_seq = feature.qualifiers["translation"][0]  # Extract protein sequence

                        # Create a unique name using start and end positions
                        feature_name = f"{start}_{end}_protein"

                        # Write to FASTA
                        fasta_out.write(f">{feature_name}\n{protein_seq}\n")
                    else:
                        print(f"Warning: No translation found for CDS at {start}-{end}")



