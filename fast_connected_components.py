# -*- coding: utf-8 -*-
import time
import os,sys
import itertools
import math
import argparse
import numpy as np
from multiprocessing import Pool
from hashlib import sha1
from datasketch import MinHash
import random, struct
from random import sample,choice
from sklearn import metrics


#############
# functions #
#############

def get_clusters(fn):
    with open(fn,'r') as f:
        f.next()#skip header
        for line in f:
            a=line.split(',')
            yield a[0],a[2]

def get_lsh(sig,nbands):
    for i,band in enumerate(np.array_split(sig,nbands)):
        yield sha1("ab" + unicode(band) + "ba"+unicode(i)).digest()
         
def get_bandwidth(n, tr):
        """
        Threshold tr = (1/b) ** (1/r) where
        b #bands
        r #rows per band
        n = b * r  #elements in signature
        """
        best = n, 1
        minerr  = float("inf")
        for r in xrange(1, n + 1):
            try:
                b = 1. / (tr ** r)
            except: 
                return best
            err = abs(n - b * r)
            if err < minerr:
                best = r
                minerr = err
        return best


def connected(seed,lshdict,doc2lsh,t):
    stk=[]
    sigs=doc2lsh[seed]
    #get lsh band signatures for seed and 
    #find match candidates
    cluster=set([seed])
    #get candidates and flatten list
    candidates=set(itertools.chain.from_iterable([lshdict[sig] for sig in doc2lsh[seed]]))
    done=candidates
    stk=candidates
    if len(stk) > 1:#candidates contain more than seed itself
        while len(stk)>0:
            cand=stk.pop()
            if cand in cluster:continue#don't check if we've already seen this 
            m1=hashcorp[cand]
            e=0
            for doc in cluster:
                m2=hashcorp[doc]
                if m2.jaccard(m1) >=t:
                    e=1
                    break
                    #we could break here if we just needed the connected component. But we need the full number of edges so we'll continue.
            if e>0:
                cluster.add(cand)
                candidates=set(itertools.chain.from_iterable([lshdict[sig] for sig in doc2lsh[cand]]))
                #update stk with candidates that haven't been seen yet. 
                stk.update(candidates.difference(done))
                done.update(candidates)
    #all candidates have been checked, full connected component is resolved. 
    return cluster 
    
def compute_clusters(obj):
    thr=obj[0]
    bandwidth=get_bandwidth(num_permutations, thr)#r
    bands=int(math.ceil(float(num_permutations)/float(bandwidth)))#b
    print "starting calculations for threshold "+str(thr)+"\nnumber of lsh bands: "+str(bands)
    sys.stdout.flush()

    start_time = time.time()
    doc_to_lsh={}
    lsh_dict={}

    for key,m in hashcorp.iteritems():
        #compute lsh 
        signatures = [sig for sig in get_lsh(m.hashvalues,bands)]
        #store signatures for this document
        doc_to_lsh[key]=signatures
        #store lsh signature to key
        for sig in signatures:
            if sig in lsh_dict:
                lsh_dict[sig].append(key)
            else:
                lsh_dict[sig]=[key]
    print("Calculating lsh signatures for threshold "+str(thr)+" took\n ---%s seconds ---\n" % (time.time() - start_time))
    sys.stdout.flush()

    #compute connected components
    start_time = time.time()
    doc2cluster={}
    cluster_info={}
    count=0

    for doc in hashcorp:
        if doc not in doc2cluster:
            cl,ed=connected(doc,lsh_dict,doc_to_lsh,thr)
            doc2cluster.update({i:count for i in cl })
            cluster_info[count]=(len(cl),ed)
            count+=1
    print("Computing connected components for threshold: "+str(thr)+" took\n--- %s seconds ---\n" % (time.time() - start_time))
        
    print "write results to file"
    start_time = time.time()
    f=open('outV3/ad2cluster_'+num_lines+'_'+str(thr)+'_'+'.csv','w')
    f.write('ad,cluster\n')
    for key, value in ad2cluster.iteritems():
        f.write(str(key)+','+str(value)+'\n')
    f.close()
    f=open('outV3/cluster_info_'+num_lines+'_'+str(thr)+'.csv','w')
    f.write('cluster,size,edges\n')
    for key, value in cluster_info.iteritems():
        f.write(str(key)+','+','.join(map(str,value))+'\n')
    f.close()
    print("Writing results to files for threshold "+str(thr)+" took:\n--- %s seconds ---\n" % (time.time() - start_time))
    
                
#Set up command line arguments
parser = argparse.ArgumentParser(description='run minhash ER experiment with given threshold')
parser.add_argument("-t", dest="threshold",type=float,help="threshold for ER", metavar="T")
parser.add_argument("-lt", dest="lt",type=float,help="lower threshold for ER", metavar="TL")
parser.add_argument("-ut", dest="ut",type=float,help="upper threshold for ER", metavar="TU")
parser.add_argument("-out", dest="out",help="output directory", metavar="OUT")
parser.add_argument("-steps", dest="steps",type=float,help="number of steps between lower and upper threshold", metavar="TSTEP")
parser.add_argument("-sigl", dest="num_permutations",type=int,help="minhash signature length", metavar="SIG")
parser.add_argument("-s", dest="suffix",help="file suffix", metavar="S")
parser.add_argument("-infile", dest="infile",help="input file",required=True, metavar="IF")
parser.add_argument('-match', dest='match', action='store_true')
parser.add_argument('-header', dest='header', action='store_true')
parser.add_argument("-numl", dest="num_lines", required=False,help="number of lines to use", metavar="NUML")
parser.add_argument("-p", dest="nump", required=False,type=int,help="number of processes for multithreading", metavar="NUMP")
parser.set_defaults(match=False)
parser.set_defaults(header=True)
parser.set_defaults(threshold=None)
parser.set_defaults(num_permutations=100)
parser.set_defaults(lt=0.0)
parser.set_defaults(ut=1.0)
parser.set_defaults(steps=2)
parser.set_defaults(nump=1)
parser.set_defaults(suffix='')
parser.set_defaults(out='out')

if __name__ == "__main__":    
    #fetch command line arguments
    args = parser.parse_args()
    num_processes=args.nump
    calc_match=args.match
    calc_clusters=args.clusters
    suffix=args.suffix
    num_lines=args.num_lines
    num_permutations=args.num_permutations#length of minhash signature

    #create output directory if it does not exist
    outdir=args.out
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    thresholds=[]
    lt=args.lt
    ut=args.ut
    steps=args.steps
    if args.threshold is not None:
        thresholds=[args.threshold]
    else:
        if None in [lt,ut,steps]: 
            print "need lower threshold, upper threshold, and number of steps"
            exit()
        else:
            thresholds=np.linspace(lt, ut, num=steps)

    #load text. Flat file for now
    print 'load text'
    start_time = time.time()
    with open(args.infile,'r') as f:
        if args.header:
            f.next()
        #TODO test robustness
        #mycorpus=[(i,set(line.encode('utf8', 'ignore').lower().split())) for i,line in enumerate(f)]
        mycorpus=[(i,set(line.lower().split())) for i,line in enumerate(f)]

    print("--- %s seconds ---" % (time.time() - start_time))

    print 'Calculate minhash signatures'
    start_time = time.time()

    #prepare dictionary of hashes
    hashcorp=dict.fromkeys([tup[0] for tup in mycorpus])
    #compute hashes
    for key,doc in mycorpus:
        #compute minhash signature
        m=MinHash(num_perm=num_permutations)
        for token in doc: m.digest(sha1(token))
        hashcorp[key]=m
    print("--- %s seconds ---" % (time.time() - start_time))

    p=Pool(num_processes)
    assignment=[ (x,) for x in thresholds]
    p.map(compute_clusters,assignment)

