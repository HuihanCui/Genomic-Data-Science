def constructEmpty(reads):
    dic = {}
    s = set()
    for read in reads:
        for i in range(71):
            kmer = read[i:i+30]
            if kmer not in s:
                s.add(kmer)
    for kmer in s:
        dic[kmer] = set()
    return dic

def constructDic(reads, dic):
    for read in reads:
        for i in range(71):
            kmer = read[i:i+30]
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
    for read in reads:
        keypart = read[-30:]
        relates = dic[keypart]
        for relateread in relates:
            if relateread != read:
                olen = overlap(read, relateread,30)
                if olen >0:
                    olaps[(read,relateread)] = olen
    return olaps




#904746
#386374


d = {'a':'t','asdsasd':'sdfsdf'}
print(len(d))