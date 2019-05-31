#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
fastq-dump --gzip --split-files SRR5930526

# De novo metagenome assembly
# SPAdes - http://cab.spbu.ru/software/spades/
spades.py -1 SRR5930526_1.fastq.gz -2 SRR5930526_2.fastq.gz -o P24.SPAdes --meta

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P24.gene P24.SPAdes/scaffolds.fasta SRR5930526_1.fastq.gz,SRR5930526_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P24.gene.region.abundance.txt P24.SPAdes/scaffolds.fasta centrifuge/data/p_compressed > P24.gene_taxon.region.abundance.txt
