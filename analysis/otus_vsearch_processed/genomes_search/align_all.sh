#!/bin/bash

#===============================================#
# Align all OTUs to HL samples
#
# March 25 2016
#===============================================#

# This alignment does not use the clustered HL.fna.
# Instead, it uses the full HL library, and reports
# multiple hits. This works well because the HL
# library actually includes multiple copies of the
# same gene from different sources. 
# 

#vsearch --usearch_global ../rep_set.fna \
#--db ../genomes/HL_16S_lib.derep.fna \
#--id 0.97 --sizeout --threads 3 --output_no_hits \
#--uc tophit_97_HL_16S_lib.uc \
#--alnout tophit_97_HL_16S_lib.aln

vsearch --usearch_global ../rep_set.fna \
--db ../genomes/HL_16S_lib.fna \
--id 1.0 --sizeout --threads 3 \
--output_no_hits --uc_allhits --maxaccepts 0 --maxrejects 4 \
--uc allhits_100_HL_16S_lib.uc \
--alnout allhits_100_HL_16S_lib.aln

vsearch --usearch_global ../rep_set.fna \
--db ../genomes/HL_16S_lib.fna \
--id 0.97 --sizeout --threads 3 \
--output_no_hits --uc_allhits --maxaccepts 0 --maxrejects 16 \
--uc allhits_97_HL_16S_lib.uc \
--alnout allhits_97_HL_16S_lib.aln
