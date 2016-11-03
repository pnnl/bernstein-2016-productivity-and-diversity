#!/bin/bash


#===============================================#
# Join all demultiplexed fastq files
#===============================================#

# Setting up folders
homedir=/people/bris469/data/hans
input=$homedir/raw_data
output=$homedir/join
rm -rf $output
mkdir $output


# Print the part of the read name that changes
cd $input
ls *f.fastq* | cut -d. -f1 > $output/readnames.txt
#ls R1*fastq | cut -c4- > $output/readnames.txt
echo Head of readnames.txt: 
head $output/readnames.txt


which vsearch

while read readname; do
	#echo pairing $readname
	vsearch --fastq_mergepairs $readname.f.fastq.gz \
		--reverse $readname.r.fastq.gz \
		--fastqout $output/$readname.m.fastq \
		--fastq_allowmergestagger --fastq_minovlen 290 &
	sleep 1
#	fastq-join -p 70 -m 290 \
#		R1.$readname R2.$readname \
#		-o $output/$readname.%.fastq \
#		-r $output/$readname.overlaps \
#		> $output/$readname.stdout &
done < $output/readnames.txt 

#cd $output
#rm -rf *.un1.*
#rm -rf *.un2.*

echo Finished!

