#!/bin/bash


#===============================================#
# Join all demultiplexed fastq files
#===============================================#

# Setting up folders
homedir=/people/bris469/data/hans
input=$homedir/join
output=$homedir/labeled
rm -rf $output
mkdir $output


# Print read names
cd $input
ls *fastq > $output/readnames.txt
echo "head of read names.txt:"
head $output/readnames.txt


which vsearch

# filter, trim, and convert to fasta...
while read readname; do
	outputname="$(echo $readname | cut -d '.' -f1)"
	vsearch --fastq_filter $input/$readname \
	--fastaout $output/$outputname.fasta \
	--fastq_maxee .5 --fastq_stripleft 24 --fastq_trunclen 253 &
	# These trimming parameters were chosen after strange reads were observed
	# during dereplication. These setting attempt capture the region 
	# between the 16S V4 primers.
	sleep 1
done < $output/readnames.txt 

#cd $output
#rm -rf *.un1.*
#rm -rf *.un2.*

sleep 5
echo Finished!

cd $output
ls *fasta > new_readnames.txt
echo New file names to add to mapping file:
cat $output/new_readnames.txt

# add the new fasta names to your mapping file under the column named "InputFileName

# ...then add labels with qiime
# validate_mapping_file.py -m map.txt -o validate_mapping -b -j InputFileName
# less validate_mapping/map.log 
# add_qiime_labels.py -m map.txt -i labeled/ -c InputFileName


