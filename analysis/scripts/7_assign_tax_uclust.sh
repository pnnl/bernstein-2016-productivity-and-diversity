#!/bin/bash

#===============================================#
# Build assign tax to a rep_set.fna using qiime's 
# parallel_assign_taxonomy_uclust.py and add
# those assignments to a .biom table
#===============================================#

# Set up directories and clean outputs
homedir=/people/bris469/data/hans
cd $homedir/otus_vsearch
rm -rf uclust_assign_tax
rm -r otu_table_mc2_w_tax*


# Make sure you have qiime loaded
#source activate qiime-git

# assign tax to the same chimera checked rep set you mapped reads back on to.
echo "running parallel_assign_taxonomy_uclust.py"
parallel_assign_taxonomy_uclust.py -i rep_set.fna \
-o uclust_assign_tax/ -T --jobs_to_start 5

# add to biome file
echo "Done!"
echo "adding metadata"
biom add-metadata -i otu_table.biom -o otu_table_w_tax.biom \
--observation-metadata-fp uclust_assign_tax/rep_set_tax_assignments.txt \
--sc-separated taxonomy --observation-header OTUID,taxonomy

biom summarize-table -i otu_table_w_tax.biom -o otu_table_w_tax.summary.txt

biom convert -i otu_table_w_tax.biom -o otu_table_w_tax.txt \
--to-tsv --header-key taxonomy


echo "Done!"
