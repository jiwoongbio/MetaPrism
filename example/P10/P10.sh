#!/bin/env bash

# Download Shotgun metagenomic sequencing data
# SRA Toolkit - https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
#fastq-dump --gzip --split-files SRR2048044
#fastq-dump --gzip --split-files SRR2048045

# De novo metagenome assembly
# SPAdes - http://cab.spbu.ru/software/spades/
#spades.py -1 SRR2048044_1.fastq.gz -2 SRR2048044_2.fastq.gz -1 SRR2048045_1.fastq.gz -2 SRR2048045_2.fastq.gz -o P10.SPAdes --meta

# Gene annotation and abundance quantification
perl ../../MetaPrism_gene.pl -p 8 P10.gene P10.SPAdes/scaffolds.fasta SRR2048044_1.fastq.gz,SRR2048044_2.fastq.gz SRR2048045_1.fastq.gz,SRR2048045_2.fastq.gz

# Taxon annotation
# Centrifuge - https://ccb.jhu.edu/software/centrifuge/
# ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed_2018_4_15.tar.gz
perl ../../MetaPrism_taxon_centrifuge.pl -p 8 P10.gene.region.abundance.txt P10.SPAdes/scaffolds.fasta centrifuge/data/p_compressed > P10.gene_taxon.region.abundance.txt

