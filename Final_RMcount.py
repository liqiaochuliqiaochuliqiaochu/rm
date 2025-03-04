from Bio import Seq


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
f1 = open("FinalRM.txt", 'w')
f1.writelines("name" + "\t" + "TypeI" + "\t" + "TypeII" + "\t" + "TypeIII" + "\t" + 'TypeIV' + "\t" + 'Total')
f1.writelines("\n")
ref = read_content("/home5/qli37/24RM/final/dna.txt")
# f2 = open("FinalRM_identity.txt",'w')
t = ''
st = ''
ty = ''
for i in range(0, len(list)):
    try:
        name = list[i][0].split(".gbff")[0]
        arr = []
        seq_name = []
        protein_result = read_content("BLASTP/" + name + "_clean.blast")
        ip = read_content("aaIn/" + name + "_clean.fasta")
    except:
        print(name)
        continue
    for j in range(0, len(protein_result)):
        try:
            # &(int(protein_result[j][3])/int(protein_result[j][5]) > 0.9)
            if ((j == 0) | (protein_result[j][0] not in seq_name)) & (float(protein_result[j][2]) > 75):
                ####Find RM system type from the DNA stock file
                t = ''
                st = ''
                for a in range(0, len(ref) - 1, 2):
                    info = ref[a][0].split(">REBASE:")[1]
                    ref_name = info.split("_")[0]
                    if protein_result[j][1].split("_")[1] == ref_name:
                        # print(info)
                        try:
                            t = info.split("Type_")[1]
                            t = t.split("_")[0]
                        except:
                            if "orphan" in info:
                                t = 'orphan'
                            else:
                                t = "Unknown"
                        # print("t,",t)
                        if "/" in ref_name:
                            st = "RM"
                        elif "." in ref_name:
                            st = ref_name.split(".")[0]
                            if 'M' in st:
                                st = 'M'
                            elif "S" in st:
                                st = "S"
                            elif "R" in st:
                                st = "R"
                            elif "C" in st:
                                st = 'C'
                                t = "III"
                        elif "restriction" in info:
                            st = "R"
                        elif "methyltransferase" in info:
                            st = "M"
                        elif "control" in info:
                            st = "C"
                            t = "III"
                            # print(protein_result[j][1])
                        else:
                            print(info)
                        break

                if (t == "Unknown") & (protein_result[j + 1][0] == protein_result[j][0]):
                    print(protein_result[j],protein_result[j+1])
                    for a in range(0, len(ref) - 1, 2):
                        info = ref[a][0].split(">REBASE:")[1]
                        ref_name = info.split("_")[0]
                        if protein_result[j + 1][1].split("_")[1] == ref_name:
                            # print(info)
                            try:
                                print(ref[a])
                                t = info.split("Type_")[1]
                                t = t.split("_")[0]
                            except:
                                if "orphan" in info:
                                    t = 'orphan'
                                else:
                                    t = "Unknown"
                            break
                        # st = "R"
                    # print("st,",st)

                arr.append([protein_result[j][0], protein_result[j][1], t, st,
                            protein_result[j][2] + "_" + protein_result[j][3] + "_" + protein_result[j][4] + "_" +
                            protein_result[j][5]])
                seq_name.append(protein_result[j][0])

        except:
            continue

    f = open("FinalRM/" + name.split(".fatsa")[0] + ".txt", 'w')
    IM = 0
    IS = 0
    IR = 0
    IRM = 0
    IIM = 0
    IIR = 0
    IIRM = 0
    IIS = 0
    IIIM = 0
    IIIS = 0
    IIIR = 0
    IIIRM = 0
    typeI = 0
    typeII = 0
    typeIII = 0
    typeIV = 0
    IVR = 0
    Other = 0
    orphan1 = 0
    orphan2 = 0
    orphan3 = 0
    orphan4 = 0
    for j in range(0, len(arr)):
        f.writelines(arr[j][0] + "\t" + arr[j][1] + "\t" + arr[j][2] + arr[j][3] + "\t" + arr[j][4])
        f.writelines("\n")
        ty = str(arr[j][2] + arr[j][3])
        try:
            if ty == "IM":
                IM = IM + 1
            elif ty == "IS":
                IS = IS + 1
            elif ty == "IR":
                IR = IR + 1
            elif ty == "IRM":
                IRM = IRM + 1
            elif ty == "IIS":
                IIS = IIS + 1
            elif ty == "IIR":
                IIR = IIR + 1
            elif ty == "IIM":
                IIM = IIM + 1
            elif "IIG" in ty:
                IIRM = IIRM + 1
            elif ty == "IIIS":
                IIIS = IIIS + 1
            elif ty == "IIIR":
                IIIR = IIIR + 1
            elif ty == "IIIM":
                IIIM = IIIM + 1
            elif ty == "IIIRM":
                IIIRM = IIIRM + 1
            elif ty == "IVR":
                IVR = IVR + 1
            # elif "III" in type:
            #     typeIII = typeIII + 1
            # elif "IV" in type:
            #     typeIV = typeIV + 1
            else:
                Other = Other + 1
        except:
            print(arr[j])
            print(name)
            print(ty)

    if IR == 0:
        orphan1 = orphan1 + IM + IS
    else:
        typeI = min(IR, IM, IS)
    if IIR == 0:
        orphan2 = orphan2 + IIM + IIS
    else:
        typeII = IIR

    if IIIR == 0:
        orphan3 = orphan3 + IIIM + IIIS
    else:
        typeIII = min(IIIM, IIIR)
        # if IIIR <= IIIM:
        #     orphan3 = orphan3 + (IIIM - IIIR)
        #     typeIII = typeIII + IIIR
        # elif IIIR > IIIM:
        #     orphan3 = orphan3 + (IIIR - IIIM)
        #     typeIII = typeIII + IIIM
        #     print(name + "!!!")

    if IVR == 0:
        orphan4 = orphan4
    else:
        typeIV = typeIV + IVR

    f1.writelines(
        name.split(".fasta")[0] + "\t"
        + str(typeI + IRM) + "\t"
        + str(typeII + IIRM) + "\t"
        + str(typeIII + IIIRM) + "\t"
        + str(typeIV) + "\t"
        + str(typeI + IRM + typeII + IIRM + typeIII + IIIRM + typeIV))
    f1.writelines("\n")
