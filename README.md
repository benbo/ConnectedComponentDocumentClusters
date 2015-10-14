# FastDocumentClusters
Compute fast document clusters for text data using MinHash.

Clusters can be computed in parallel for different Jaccard similarity thresholds. The parallel implementation however is naive and potentially requires lots of memory.  

##Dependencies
numpy, scikit-learn, datasketch

to install datasketch (https://github.com/ekzhu/datasketch) do

    pip2.7 install datasketch -U


##Example
    python2.7 fast_document_clusters.py -infile test.txt -lt 0.9 -ut 1.0 -steps 2 -p 2 -suff test
