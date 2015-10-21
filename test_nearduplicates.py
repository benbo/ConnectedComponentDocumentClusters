import nearduplicates
from hashlib import sha1
from datasketch import MinHash

corpus=[
        'Vice President Biden announced that he won\'t enter the race for the 2016 presidential nomination',
        'Vice President Joe Biden said that he will no enter the race for the 2016 presidential nomination',
        'When 2012 GOP nominee Mitt Romney picked Ryan to be his vice presidential pick, the choice was largely seen as a risky one',
        'When 2012 republican nominee Romney picked Ryan to be his vice presidential pick, it was largely seen as a risky choice'
        ]

def test():
    hashcorp=nearduplicates.run_getminhash({'objects':[(i,x) for i,x in enumerate(corpus)]})

    print "actual_jaccard "
    for i in xrange(4):
        for j in xrange(i+1,4):
            print i,j
            s1 = set(corpus[i])
            s2 = set(corpus[j])
            print "jaccard    : " + str( float(len(s1.intersection(s2)))/float(len(s1.union(s2))))
            a1= hashcorp[i]
            a2= hashcorp[j]
            print "minhash jac: " + str(nearduplicates.run_jaccard_array({'signatures':(a1,a2)}))
            m1, m2 = MinHash(), MinHash()
            for d in corpus[i].split():
                m1.digest(sha1(d.encode('utf8')))
            for d in corpus[j].split():
                m2.digest(sha1(d.encode('utf8')))
            print("Estimated Jaccard for data1 and data2 is", m1.jaccard(m2))



if __name__ == "__main__":
    test()
