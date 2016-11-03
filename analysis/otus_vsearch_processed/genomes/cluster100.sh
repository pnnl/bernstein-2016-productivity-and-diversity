#!/bin/bash

#===============================================#
# Cluster reads at 100% (substring) 
#
# March 25 2016
#===============================================#

vsearch --cluster_fast HL_16S_lib.fna \
--id 1 --sizeout --threads 2 \
--centroids HL_16S_lib.cluster100.fna --uc HL_16S_lib.cluster100.uc