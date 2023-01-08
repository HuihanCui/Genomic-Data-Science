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
        for i in range(len(read)-10+1):
            kmer = read[i:i+10]
            if kmer not in s:
                s.add(kmer)
    for kmer in s:
        dic[kmer] = set()
    return dic,list(s)

def constructDic(reads, dic):
    for read in reads:
        for i in range(len(read)-10+1):
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

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    lists = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            lists = []
        if len(sup) == len(shortest_sup):
            lists.append(sup)
    return shortest_sup,lists  # return shortest

def pick(reads,results):
    reada, readb = None, None
    best_olen = 0
    for (a,b) in results.keys():
        olen = results[(a,b)]
        if olen > best_olen:
            reada, readb = a,b 
            best_olen = olen
    return reada, readb, best_olen

def greedy(reads):
    empty,lists = constructEmpty(reads)
    real = constructDic(reads,empty)
    result,count = betteroverlap(reads,real)
    
    reada, readb, olen = pick(reads,result)
    while olen>0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[olen:])
        
        empty,lists = constructEmpty(reads)
        real = constructDic(reads,empty)
        result,count = betteroverlap(reads,real)
        
        reada, readb, olen = pick(reads,result)
    return ''.join(reads)

#strs = ['CCT', 'CTT', 'TGC','TGG','GAT','ATT']
#a,b = scs(strs)
#print(a)
#print(len(b))

seqs, _ = readFastq('ads1_week4_reads.fq')


a = greedy(seqs)

print(len(a))
print(a.count('A'))
print(a.count('T'))

