#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
fastq-dump --gzip --split-files --log-level err SRR5930516

# De novo metagenome assembly
# MEGAHIT - https://github.com/voutcn/megahit
megahit -1 SRR5930516_1.fastq.gz -2 SRR5930516_2.fastq.gz -o P52.megahit

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P52.gene P52.megahit/final.contigs.fa SRR5930516_1.fastq.gz,SRR5930516_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P52.gene.region.abundance.txt P52.megahit/final.contigs.fa centrifuge/data/p_compressed > P52.gene_taxon.region.abundance.txt

