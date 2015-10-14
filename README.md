# FastConnectedComponents
Compute fast connected components for text data using MinHash.

Connected components can be computed in parallel for different Jaccard similarity thresholds. The parallel implementation however is naive and potentially requires lots of memory.  

##Dependencies
numpy, scikit-learn, datasketch

to install datasketch (https://github.com/ekzhu/datasketch) do

    pip2.7 install datasketch -U


##Example
    python2.7 fast_connected_components.py -infile test.txt -lt 0.9 -ut 1.0 -steps 2 -p 2 -suff test
