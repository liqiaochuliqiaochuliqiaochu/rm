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


list = read_content("3.txt")
f = open("BLASTP.sh", "w")
f.writelines("#!/bin/bash")
f.writelines("\n")
for i in range(0, len(list)):
    name = list[i][0].split(".gbff")[0]
    f.writelines("blastp -query aaIn/"+name +"_clean.fasta -out BLASTP/" + name +"_clean.blast -db aaIn/aa.blastdb -outfmt \"6 qseqid sseqid pident length qlen slen evalue bitscore\" -evalue 1e-30 -num_threads 12 ")
    f.writelines("\n")
