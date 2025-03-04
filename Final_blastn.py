def read_content(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    i = 0
    data = []
    for line in lines:
        line = line.strip().split('\t')  # strip \n, split by space
        data.append(line)  # record data
        i = i + 1
    return data


list = read_content("list.txt")
f = open("BLASTN.sh", "w")
f.writelines("#!/bin/bash")
f.writelines("\n")
for i in range(0, len(list)):
    name = list[i][0].split(".gbff")[0]
    f.writelines(
        "blastn -query In/" + name + ".fasta -out BLASTN/" + name + ".blast -db /home5/qli37/24RM/final/test.blastdb -outfmt \"6 qseqid sseqid pident length qlen slen evalue bitscore\" -evalue 1e-30 -num_threads 12 -qcov_hsp_perc 99 -perc_identity 95 -max_target_seqs 50")
    f.writelines("\n")
