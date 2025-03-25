def read_content(filename):
    """
    Reads a tab-delimited file and returns its contents as a list of lists.
    
    :param filename: Path to the file to be read.
    :return: A list of lists, where each sublist represents a line in the file.
    """
    with open(filename, 'r') as file:
        return [line.strip().split('\t') for line in file]

def generate_cd_hit_script(input_file, output_script):
    """
    Generates a shell script to run CD-HIT for sequence clustering.

    :param input_file: Path to the input file containing genome names.
    :param output_script: Path to the output shell script.
    """
    list_data = read_content(input_file)
    
    with open(output_script, "w") as f:
        f.write("#!/bin/bash\n\n")  # Shebang for bash script

        for entry in list_data:
            name = entry[0].split(".gbk")[0]  # Extract base name without .gbff extension
            
            f.write(
                f"cd-hit -i aaIn/{name}.fasta -o aaIn/{name}_clean.fasta -c 0.95\n"
            )

# Usage
generate_cd_hit_script("list.txt", "CD.sh")
