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
f = open("CD.sh", "w")
f.writelines("#!/bin/bash")
f.writelines("\n")
for i in range(0, len(list)):
    name = list[i][0].split(".gbff")[0]
    f.writelines("cd-hit -i aaIn/"+name +".fasta -o aaIn/" + name +"_clean.fasta -c 0.95")
    f.writelines("\n")
