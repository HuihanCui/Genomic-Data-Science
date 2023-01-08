def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        count = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                count += 1
            if count > 2:
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def findGCByPos(reads):
    gc = [0]*100
    total = [0]*100
    for read in reads: 
        for i in range(len(read)):
            if read[i] == 'G' or read[i] == 'C':
                gc[i] += 1
            total[i] += 1
            
    for i in range(len(gc)):
        if total[i] >0:
            gc[i] = gc[i] / float(total[i])
            
    return gc

def findBadQua(reads):
    qua = [0]*100
    for read in reads:
        for i in range(len(read)):
            qua[i] += ord(read[i])-33
    return qua

#genome = readGenome('lambda_virus.fa')
seqs, quas = readFastq('ERR037900_1.first1000.fastq')

scores = findGCByPos(seqs)
for i in range(len(scores)):
    if scores[i] == 0.596:
        pass
    
qua_s = findBadQua(quas)
for i in range(len(qua_s)):
    if qua_s[i] == 4526:
        pass

print(qua_s[60:80])
    
#print(naive(reverseComplement('AGTCGA'),genome))