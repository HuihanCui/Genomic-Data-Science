file = open("dna2.fasta")

def num_records():
    count = 0
    for line in file:
        line = line.rstrip()
        if line[0] == ">":
            count+=1
    return count
        
def lengths():
    lengths = []
    seqs = {}
    for line in file:
        line = line.rstrip()
        if line[0] == ">":
            words = line.split()
            name = words[0]
            seqs[name] = ''
        else:
            seqs[name] = seqs[name] + line
    for i in seqs.values():
        lengths.append(len(i))
    return lengths
        
def seqs():
    seqs = {}
    file.seek(0)
    for line in file:
        line = line.rstrip()
        if line[0] == ">":
            words = line.split()
            name = words[0]
            seqs[name] = ''
        else:
            seqs[name] = seqs[name] + line
    return seqs


def orf(list):
    result = []
    dic = {}
    start = False
    for i in range(2,len(list),3):
        codon = list[i:i+3].lower()
        if start == False:
            if codon == "atg":
                dic["start"]= i
                start = True
        else:
            if codon == "taa" or codon == "tag" or codon == "tga":
                dic["end"] = i
                dic["length"] = dic["end"] - dic["start"]
                result.append(dic["length"])
                start = False
                dic = {}
    if len(result) == 0:
        return 0
    return max(result)

hhh =  list(seqs().values())

def repeats(str):
    sum = 0
    for i in hhh:
        sum = sum + i.count(str)
    return sum
 
orfs = []           
       
for i in range(len(hhh)):
     
     orfs.append(orf(hhh[i]))
    
answer = list(seqs().keys())
for i in range(len(answer)):
    if answer[i] == ">gi|142022655|gb|EQ086233.1|16":
        #print(i)
        pass
       
       
def all12():
    choice = ['A','T','G','C']
    result = []
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    for e in range(4):
                        for f in range(4):
                            for g in range(4):
                                for h in range(4):
                                    for i in range(4):
                                        for j in range(4):
                                            for k in range(4):
                                                for l in range(4):
                                                    result.append(choice[a]+choice[b]+choice[c]+choice[d]+choice[e]+choice[f]+choice[g]+choice[h]+choice[i]+choice[j]+choice[k]+choice[l])
    return result

repeatss = []
for i in all12():
    repeatss.append(repeats(i))

print(max(repeatss))     
    