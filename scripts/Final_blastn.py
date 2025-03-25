def read_content(filename):
    """
    Reads a tab-delimited file and returns its contents as a list of lists.
    
    :param filename: Path to the file to be read.
    :return: A list of lists, where each sublist represents a line in the file.
    """
    with open(filename, 'r') as file:
        return [line.strip().split('\t') for line in file]

def generate_blast_script(list_file, output_script, db_name="rebase_dna"):
    """
    Generates a shell script to run BLASTN on multiple query files.

    :param list_file: Path to the input file containing list of genome names.
    :param output_script: Path to the output shell script file.
    :param db_name: Name of the BLAST database.
    """
    list_data = read_content(list_file)
    
    with open(output_script, "w") as f:
        f.write("#!/bin/bash\n\n")  # Shebang for bash script

        for entry in list_data:
            name = entry[0].split(".gbk")[0]  # Extract base name without .gbff extension
            
            f.write(
                f"blastn -query In/{name}.fasta "
                f"-out BLASTN/{name}.blast "
                f"-db {db_name} "
                f"-outfmt \"6 qseqid sseqid pident length qlen slen evalue bitscore\" "
                f"-evalue 1e-30 -num_threads 12 -qcov_hsp_perc 99 -perc_identity 95 "
                f"-max_target_seqs 50\n"
            )

# Usage
generate_blast_script("list.txt", "BLASTN.sh")
