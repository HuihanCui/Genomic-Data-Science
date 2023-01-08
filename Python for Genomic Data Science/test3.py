file = open("dna2.fasta")

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

hhh =  list(seqs().values())

for i in hhh:
    bool = i.isupper()
    print(bool)