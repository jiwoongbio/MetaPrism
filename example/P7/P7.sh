#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
#fastq-dump --gzip --split-files SRR2047840
#fastq-dump --gzip --split-files SRR2047841

# De novo metagenome assembly
# SPAdes - http://cab.spbu.ru/software/spades/
#spades.py -1 SRR2047840_1.fastq.gz -2 SRR2047840_2.fastq.gz -1 SRR2047841_1.fastq.gz -2 SRR2047841_2.fastq.gz -o P7.SPAdes --meta

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P7.gene P7.SPAdes/scaffolds.fasta SRR2047840_1.fastq.gz,SRR2047840_2.fastq.gz SRR2047841_1.fastq.gz,SRR2047841_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P7.gene.region.abundance.txt P7.SPAdes/scaffolds.fasta centrifuge/data/p_compressed > P7.gene_taxon.region.abundance.txt

