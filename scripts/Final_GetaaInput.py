from Bio import SeqIO
from Bio.Seq import Seq


def read_content(filename):
    """
    Reads a tab-delimited file and returns its contents as a list of lists.
    Each line is stripped of newline characters and split into a list by tab.
    
    :param filename: Path to the file to be read.
    :return: A list of lists, where each sublist represents a line in the file.
    """
    file = open(filename, 'r')
    lines = file.readlines()
    data = []
    
    for line in lines:
        line = line.strip().split('\t')  # Strip newline and split by tab
        data.append(line)  # Store the processed line
    
    return data

list = read_content("list.txt")
ref = read_content("protein_seqs.fasta")
ref_name = ""
ref_seq = ''
tempname = ''
name_arr = []
protein_sequence = ''
f2 = open("aaIn/aadb.fasta", 'w')
for i in range(0, len(list)):
    name = list[i][0].split(".gbk")[0]
    f1 = open("aaIn/" + name + ".fasta", 'w')
    data = read_content("BLASTN/" + name + ".blast")
    protein_file = read_content("In/" + name + "_protein.fasta")
    fasta_file = read_content('In/' + name + ".fasta")
    in_arr = []
    for j in range(0, len(data)):
        ref_name = data[j][1].split("_")[1]
        if (data[j][0] not in in_arr) & (int(data[j][3])/int(data[j][5]) > 0.9):
            protein_sequence = ''
            for x in range(0, len(protein_file), 2):
                # Loop through features to find the matching protein sequence
                if data[j][0] in protein_file[x][0]:
                    protein_sequence = protein_file[x + 1][0]
                    break
            if protein_sequence == '':
                for x in range(0, len(fasta_file), 2):
                    if data[j][0] in fasta_file[x][0]:
                        dna_sequence = fasta_file[x + 1][0]
                        dna_sequence = Seq(dna_sequence)
                        protein_sequence = str(dna_sequence.translate(to_stop=False))
                        protein_sequence = protein_sequence.replace("*", "X")
                        break
            f1.writelines(">" + data[j][0])
            f1.writelines("\n")
            f1.writelines(protein_sequence)
            f1.writelines("\n")
            in_arr.append(data[j][0])

        if (int(data[j][3])/int(data[j][5]) > 0.9) & (ref_name not in name_arr):
            for a in range(0, len(ref), 2):
                if ref_name == ref[a][0].split("_")[1]:
                    ref_seq = ref[a + 1][0]
                    name_arr.append(ref_name)
                    f2.writelines(">" + data[j][1])
                    f2.writelines("\n")
                    f2.writelines(ref_seq)
                    f2.writelines("\n")
                    break
        else:
            continue
