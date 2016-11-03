#!/bin/bash

#===============================================#
# Pull all the OTUs and their taxonomies 
# from these families.
# 
#c("o__Spirochaetales f__Spirochaetaceae", 
# "o__Oscillatoriales f__Phormidiaceae",
# "o__Pseudanabaenales f__Pseudanabaenaceae",
# "o__SBR1031 f__A4b",
# "o__Rhodospirillales",
# "o__Rhodobacterales f__Rhodobacteraceae",
# "o__Rhizobiales",
# "o__Rhodospirillales f__Rhodospirillaceae",
# "o__Bacteroidales f__ML635J-40",
# "o__Phycisphaerales")
# 
# July 10 2015
#===============================================#

# Set up directories and clean outputs
homedir=/Volumes/Projects/hans/colin/otus_vsearch

cd $homedir
rm -rf topOTUs
mkdir topOTUs

#grep "o__Spirochaetales; f__Spirochaetaceae" < otu_table_w_tax.txt | wc -l
#grep "o__Oscillatoriales; f__Phormidiaceae" < otu_table_w_tax.txt | wc -l
#grep "o__Pseudanabaenales; f__Pseudanabaenaceae" < otu_table_w_tax.txt | wc -l
#grep "o__SBR1031; f__A4b" < otu_table_w_tax.txt | wc -l
#grep "o__Rhodospirillales; f__" < otu_table_w_tax.txt | wc -l
#grep "o__Rhodobacterales; f__Rhodobacteraceae" < otu_table_w_tax.txt | wc -l
#grep "o__Rhizobiales; f__" < otu_table_w_tax.txt | wc -l
#grep "o__Rhodospirillales; f__Rhodospirillaceae" < otu_table_w_tax.txt | wc -l
#grep "o__Bacteroidales; f__ML635J-40" < otu_table_w_tax.txt | wc -l
#grep "o__Phycisphaerales; f__" < otu_table_w_tax.txt | wc -l



echo "o__Spirochaetales; f__Spirochaetaceae" >> topOTUs/OTUsintopFamilies.txt
grep "o__Spirochaetales; f__Spirochaetaceae" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Oscillatoriales; f__Phormidiaceae" >> topOTUs/OTUsintopFamilies.txt
grep "o__Oscillatoriales; f__Phormidiaceae" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Pseudanabaenales; f__Pseudanabaenaceae" >> topOTUs/OTUsintopFamilies.txt
grep "o__Pseudanabaenales; f__Pseudanabaenaceae" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__SBR1031; f__A4b" >> topOTUs/OTUsintopFamilies.txt
grep "o__SBR1031; f__A4b" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Rhodospirillales; f__" >> topOTUs/OTUsintopFamilies.txt
grep "o__Rhodospirillales; f__;" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Rhodobacterales; f__Rhodobacteraceae" >> topOTUs/OTUsintopFamilies.txt
grep "o__Rhodobacterales; f__Rhodobacteraceae" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Rhizobiales; f__" >> topOTUs/OTUsintopFamilies.txt
grep "o__Rhizobiales; f__;" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Rhodospirillales; f__Rhodospirillaceae" >> topOTUs/OTUsintopFamilies.txt
grep "o__Rhodospirillales; f__Rhodospirillaceae" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Bacteroidales; f__ML635J-40" >> topOTUs/OTUsintopFamilies.txt
grep "o__Bacteroidales; f__ML635J-40" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt

echo "o__Phycisphaerales; f__" >> topOTUs/OTUsintopFamilies.txt
grep "o__Phycisphaerales; f__;" < otu_table_w_tax.txt >> topOTUs/OTUsintopFamilies.txt





