# FastConnectedComponents
Compute fast connected components for text data using MinHash.

Connected components can be computed in parallel for different Jaccard similarity thresholds. The parallel implementation however is naive and potentially requires lots of memory.  

##Dependencies

https://github.com/ekzhu/datasketch
    pip2.7 install datasketch -U
