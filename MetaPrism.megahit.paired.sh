# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

codeDir=`dirname $0`
fastqFile1=$1
fastqFile2=$2
threads=$3
prefix=$4

if [ -z "$prefix" ]; then
	echo 'Usage: ./MetaPrism.megahit.paired.sh <read.1.fastq> <read.2.fastq> <threads> <prefix>' 1>&2
	exit 1
fi

# Assembly
time megahit -1 $fastqFile1 -2 $fastqFile2 -t $threads -o $prefix.megahit

# Mapping
time bwa index $prefix.megahit/final.contigs.fa
time bwa mem -t $threads -T 19 -Y $prefix.megahit/final.contigs.fa $fastqFile1 $fastqFile2 | gzip > $prefix.megahit.sam.gz

# Sorting
time samtools sort --threads $threads $prefix.megahit.sam.gz > $prefix.megahit.sorted.bam && rm $prefix.megahit.sam.gz
time samtools index $prefix.megahit.sorted.bam

# Functional and taxonomic profiling
time perl $codeDir/MetaPrism.pl -p $threads $prefix.megahit/final.contigs.fa > $prefix.megahit.MetaPrism.txt
time perl $codeDir/MetaPrism.abundance.pl -p $threads $prefix.megahit.MetaPrism.txt $prefix.megahit.sorted.bam > $prefix.megahit.MetaPrism.abundance.txt
