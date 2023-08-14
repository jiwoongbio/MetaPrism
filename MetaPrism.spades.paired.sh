# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

codeDir=`dirname $0`
fastqFile1=$1
fastqFile2=$2
threads=$3
prefix=$4

if [ -z "$prefix" ]; then
	echo 'Usage: ./MetaPrism.spades.paired.sh <read.1.fastq> <read.2.fastq> <threads> <prefix>' 1>&2
	exit 1
fi

# Assembly
time spades.py -1 $fastqFile1 -2 $fastqFile2 -o $prefix.spades --threads $threads

# Mapping
time bwa index $prefix.spades/scaffolds.fasta
time bwa mem -t $threads -T 19 -Y $prefix.spades/scaffolds.fasta $fastqFile1 $fastqFile2 | gzip > $prefix.spades.sam.gz

# Sorting
time samtools sort --threads $threads $prefix.spades.sam.gz > $prefix.spades.sorted.bam && rm $prefix.spades.sam.gz
time samtools index $prefix.spades.sorted.bam

# Functional and taxonomic profiling
time perl $codeDir/MetaPrism.pl -p $threads $prefix.spades/scaffolds.fasta > $prefix.spades.MetaPrism.txt
time perl $codeDir/MetaPrism.abundance.pl -p $threads $prefix.spades.MetaPrism.txt $prefix.spades.sorted.bam > $prefix.spades.MetaPrism.abundance.txt
