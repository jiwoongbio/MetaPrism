#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
fastq-dump --gzip --split-files SRR5930522

# De novo metagenome assembly
# SPAdes - http://cab.spbu.ru/software/spades/
spades.py -1 SRR5930522_1.fastq.gz -2 SRR5930522_2.fastq.gz -o P32.SPAdes --meta

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P32.gene P32.SPAdes/scaffolds.fasta SRR5930522_1.fastq.gz,SRR5930522_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P32.gene.region.abundance.txt P32.SPAdes/scaffolds.fasta centrifuge/data/p_compressed > P32.gene_taxon.region.abundance.txt

