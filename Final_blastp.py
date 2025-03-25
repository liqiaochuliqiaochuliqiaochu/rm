def read_content(filename):
    """
    Reads a tab-delimited file and returns its contents as a list of lists.
    
    :param filename: Path to the file to be read.
    :return: A list of lists, where each sublist represents a line in the file.
    """
    with open(filename, 'r') as file:
        return [line.strip().split('\t') for line in file]

def generate_blastp_script(input_file, output_script, db_name="aaIn/aa.blastdb"):
    """
    Generates a shell script to run BLASTP on multiple query files.

    :param input_file: Path to the input file containing genome names.
    :param output_script: Path to the output shell script.
    :param db_name: Path to the BLAST database.
    """
    list_data = read_content(input_file)
    
    with open(output_script, "w") as f:
        f.write("#!/bin/bash\n\n")  # Shebang for bash script

        for entry in list_data:
            name = entry[0].split(".gbk")[0]  # Extract base name without .gbff extension
            
            f.write(
                f"blastp -query aaIn/{name}_clean.fasta "
                f"-out BLASTP/{name}_clean.blast "
                f"-db {db_name} "
                f"-outfmt \"6 qseqid sseqid pident length qlen slen evalue bitscore\" "
                f"-evalue 1e-30 -num_threads 12\n"
            )

# Usage
generate_blastp_script("list.txt", "BLASTP.sh")
