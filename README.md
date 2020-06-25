# cactus-connectivity
Determines how well-connected a .hal cactus graph file is by measuring how many bases are involved in alignments between each assembly in the graph.

The overall goal of cactus-connectivity is to see which bases in a given 
assembly/reference aren't aligned to other sequences.

Liftovers help with that by outputting a "target" bedfile that shows the regions in 
the target sequence that are aligned to the source sequence in positions specified by 
the source bedfile. In our case, the source bedfile simply includes _all_ positions in 
the source. Thus, we can find all the regions in the target sequence aligned to the 
source seq in the graph.

When we want to know the bases that are "free" of alignments in a certain assembly, we
want that assembly to be the target. The output of calculate_bases_unmapped will tell us
how many bases in the target aren't aligned to the source.