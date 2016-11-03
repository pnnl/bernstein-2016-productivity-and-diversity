#!/bin/bash

#===============================================#
# Using vsearch, pick OTUs
#===============================================#

# Setting up folders
homedir=/people/bris469/data/hans
output=$homedir/otus_vsearch


cd $output
rm -f rep_set.*


# clustering (uclust)
vsearch --cluster_smallmem seqs.checked_denovo.fna --threads 20 \
--id 0.97 --centroids rep_set.fna --log rep_set.log \
--sizein --xsize --usersort --relabel OTU_
#Clusters: 3877 Size min 2, max 6877690, avg 86.1

