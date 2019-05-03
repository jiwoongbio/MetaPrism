#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
#fastq-dump --gzip --split-files SRR2047825
#fastq-dump --gzip --split-files SRR2047828

# De novo metagenome assembly
# SPAdes - http://cab.spbu.ru/software/spades/
#spades.py -1 SRR2047825_1.fastq.gz -2 SRR2047825_2.fastq.gz -1 SRR2047828_1.fastq.gz -2 SRR2047828_2.fastq.gz -o P4.SPAdes --meta

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P4.gene P4.SPAdes/scaffolds.fasta SRR2047825_1.fastq.gz,SRR2047825_2.fastq.gz SRR2047828_1.fastq.gz,SRR2047828_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P4.gene.region.abundance.txt P4.SPAdes/scaffolds.fasta centrifuge/data/p_compressed > P4.gene_taxon.region.abundance.txt

