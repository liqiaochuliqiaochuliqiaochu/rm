# This script processes dna_seqs.txt and converts it into a FASTA format
# The sequence name is generated as 'order + name'

def process_sequences(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        order = 1  # Initialize order starting from 1
        lines = infile.readlines()
        sequence_name = None
        sequence = []
        
        # Skip the initial unwanted lines (if any)
        for line in lines:
            line = line.strip()
            if line.startswith("REBASE"):
                continue  # Skip the lines with 'REBASE' or any unwanted header
            
            # Look for the header line (which starts with '>')
            if line.startswith(">"):
                # If there's already a sequence, write the previous one to the file
                if sequence:
                    outfile.write(f">{sequence_name}\n")
                    outfile.write("".join(sequence) + "\n")  # Write sequence in one line
                # Extract the sequence identifier part (after '>'), skip the 'REBASE:' part
                parts = line.split("\t")
                sequence_name = f"{order}_{parts[0].split(':')[1]}"  # Using order + M.AaaS1ORF662P
                sequence = []
                order += 1  # Increment order
            else:
                # Collect the sequence (remove any whitespace)
                cleaned_line = line.replace(" ", "")  # Removing spaces from sequence
                if cleaned_line.endswith("<>"):
                    cleaned_line = cleaned_line[:-2]  # Remove the "<>" at the end
                sequence.append(cleaned_line)
                
        # Write the last sequence to the file
        if sequence:
            outfile.write(f">{sequence_name}\n")
            outfile.write("".join(sequence) + "\n")  # Write sequence in one line

# Input and output file paths
input_file = "dna_seqs.txt"
output_file = "dna_seqs.fasta"

process_sequences(input_file, output_file)

print(f"FASTA file '{output_file}' has been generated successfully.")

