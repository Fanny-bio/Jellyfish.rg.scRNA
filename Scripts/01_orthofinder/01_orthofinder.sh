ulimit -n 2600
##create a directory that contains the longest-transcripts of taxa 

##conda activate orthofinder
orthofinder -f `pwd`/longest-gene_15species -t 60 -a 60 -M msa -S diamond
