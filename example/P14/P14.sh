#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
fastq-dump --gzip --split-files --log-level err SRR5930499
fastq-dump --gzip --split-files --log-level err SRR5930500

# De novo metagenome assembly
# MEGAHIT - https://github.com/voutcn/megahit
megahit -1 SRR5930499_1.fastq.gz,SRR5930500_1.fastq.gz -2 SRR5930499_2.fastq.gz,SRR5930500_2.fastq.gz -o P14.megahit

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P14.gene P14.megahit/final.contigs.fa SRR5930499_1.fastq.gz,SRR5930499_2.fastq.gz SRR5930500_1.fastq.gz,SRR5930500_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P14.gene.region.abundance.txt P14.megahit/final.contigs.fa centrifuge/data/p_compressed > P14.gene_taxon.region.abundance.txt

