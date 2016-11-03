#!/bin/bash

#===============================================#
# Identify how many reads in our samples closely
# match our standards. This let's us estimate 
# the cross contamination rate between samples
# on this run.
#===============================================#

homedir=/people/bris469/data/hans
# Used subsampled reads
reads=$homedir/one_pct.fna
ref=$homedir/standards/Mock_Hot_Lake_Community.aligned.v4.fna
out=$homedir/standards/map
which vsearch

cd $homedir

rm -rf $out
mkdir $out

# Subsample for speed
# I did this once manually,
#vsearch -fastx_subsample combined_seqs.fna.gz -sample_pct 1 -randseed 1 -fastaout one_pct.fna 


vsearch -search_exact $reads -db $ref \
-strand plus -id 1.0 -uc $out/map_100.uc -threads 22
#Matching query sequences: 59708 of 273806 (21.81%)

vsearch -usearch_global $reads -db $ref \
-strand plus -id 0.99 -uc $out/map_99.uc -threads 22
#Matching query sequences: 87993 of 273806 (32.14%)


# Lowered then normal (97% cutoff)
vsearch -usearch_global $reads -db $ref \
-strand plus -id 0.90 -uc $out/map_90.uc -threads 22
#Matching query sequences: 164510 of 273806 (60.08%)

vsearch -usearch_global $reads -db $ref \
-strand plus -id 0.80 -uc $out/map_80.uc -threads 22
#Matching query sequences: 239413 of 273806 (87.44%)



cd $out
# make a table of OTUs
uc2otutab_mod2.py map_100.uc > map_100.txt
uc2otutab_mod2.py map_99.uc > map_99.txt
uc2otutab_mod2.py map_90.uc > map_90.txt
uc2otutab_mod2.py map_80.uc > map_80.txt




