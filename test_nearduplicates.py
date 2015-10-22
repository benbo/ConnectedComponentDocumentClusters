import nearduplicates
from hashlib import sha1

corpus=[
        'Vice President Biden announced that he will not enter the race for the 2016 presidential nomination',
        'Vice President Joe Biden said that he will not enter the race for the 2016 presidential nomination',
        'This is exactly the probability of collision we would expect if the hash function assigned truly random hash codes to every key',
        'This is exactly the probability of collision one would expect if our hash function assigned a truly random hash to every key'
        ]

def test_approximation():
    hashcorpus = [nearduplicates.run_getminhash({'id':i,'text':x}) for i,x in enumerate(corpus)]
    #create a dictionary from the nodes
    hashdict={obj['id']:obj['hashv'] for obj in hashcorpus}
    for i in xrange(4):
        for j in xrange(i+1,4):
            print i,j
            s1 = set(corpus[i].split())
            s2 = set(corpus[j].split())
            print "jaccard    : " + str( float(len(s1.intersection(s2)))/float(len(s1.union(s2))))
            a1= hashdict[i]
            a2= hashdict[j]
            print "minhash jac: " + str(nearduplicates.run_jaccard_array({'signatures':(a1,a2)}))

def test():
    #create minhash signature for all documents
    hashcorpus = [nearduplicates.run_getminhash({'id':i,'text':x}) for i,x in enumerate(corpus)]
    #get lsh buckets for threshold 0.7.
    #when searching for near duplicates, always make sure that the lsh signatures were created for the same threshold
    doc_to_lsh,lsh_dict = nearduplicates.run_lsh_batch({'threshold':0.7,'data':hashcorpus})

    #create a dictionary from the nodes so we can look up the minhash for each document id
    hashdict={obj['id']:obj['hashv'] for obj in hashcorpus}
    #find near duplicates for each text in the corpus
    for i in xrange(4):
        print 'near duplicates for: \t'+corpus[i]
        #Find near duplicates. 
        #The returned cluster is a set of ids.
        #The returned cluster currently also contains the seed itself.
        cluster = nearduplicates.run_near_duplicates({'seed':i,'hashcorp':hashdict,'doc_to_lsh':doc_to_lsh,'lsh_dict':lsh_dict,'threshold':0.7})
        for j in cluster:
            if j!=i:
                print '\t'+corpus[j]




if __name__ == "__main__":
    print "testing minhash's jaccard appoximation for short documents"
    test_approximation()
    print "test workflow for finding near duplicates"
    test()
