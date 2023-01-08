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

def constructEmpty(reads):
    dic = {}
    s = set()
    for read in reads:
        for i in range(91):
            kmer = read[i:i+10]
            if kmer not in s:
                s.add(kmer)
    for kmer in s:
        dic[kmer] = set()
    return dic,list(s)

def constructDic(reads, dic):
    for read in reads:
        for i in range(91):
            kmer = read[i:i+10]
            dic[kmer].add(read)
    return dic


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
   
def betteroverlap(reads,dic):
    olaps = {}
    count = 0
   
    for read in reads:
        find = False
        keypart = read[-10:]
        relates = dic[keypart]
        for relateread in relates:
            if relateread != read:
                olen = overlap(read, relateread,10)
                if olen >0:
                    olaps[(read,relateread)] = olen
                    find = True
        if find == True:
            count += 1
    return olaps,count
   
seqs, _ = readFastq('ads1_week4_reads.fq')
empty,lists = constructEmpty(seqs)
real = constructDic(seqs,empty)
result,count = betteroverlap(seqs,real)

results = result.keys()
for (a,b) in result:
    if (a in seqs) == False:
        print(1)
#reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
#empty,lists = constructEmpty(reads)
#real = constructDic(reads,empty)
#result = betteroverlap(reads,real)
#print(len(list(result.keys())))
#print(list(result.keys())[:10])
#print(list(result.keys())[:1])
#print(min(list(result.values())))


print('sdfsdfs'.count('s'))
