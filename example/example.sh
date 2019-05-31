#!/bin/env bash

# Prepare database files
perl ../MetaPrism_gene_prepare.pl

# Gene-taxon annotation and abundance quantification
(cd P10; ./P10.sh)
(cd P14; ./P14.sh)
(cd P23; ./P23.sh)
(cd P25; ./P25.sh)
(cd P34; ./P34.sh)
(cd P39; ./P39.sh)
(cd P8; ./P8.sh)
(cd P16; ./P16.sh)
(cd P22; ./P22.sh)
(cd P24; ./P24.sh)
(cd P30; ./P30.sh)
(cd P32; ./P32.sh)
(cd P42; ./P42.sh)

# Compare sample groups and identify differentially-abundant genes
perl ../MetaPrism_comparison.pl -F gene sample.group.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P39=P39/P39.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P16=P16/P16.gene_taxon.region.abundance.txt \
	P22=P22/P22.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P42=P42/P42.gene_taxon.region.abundance.txt \
	> gene.comparison.txt

awk -F'\t' '(NR == 1 || ($4 > 1 && $5 < 0.01))' gene.comparison.txt > gene.comparison.filtered.txt

# Generate heatmap
perl ../MetaPrism_heatmap.pl -F gene -s -g gene.comparison.filtered.txt -r both \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P39=P39/P39.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P16=P16/P16.gene_taxon.region.abundance.txt \
	P22=P22/P22.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P42=P42/P42.gene_taxon.region.abundance.txt \
	> gene.heatmap.html

# Generate table
perl ../MetaPrism_table.pl -F taxon_average -s \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P39=P39/P39.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P16=P16/P16.gene_taxon.region.abundance.txt \
	P22=P22/P22.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P42=P42/P42.gene_taxon.region.abundance.txt \
	> taxon.table.txt

# Prediction model accuracy
perl ../MetaPrism_prediction.pl sample.group.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P39=P39/P39.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P16=P16/P16.gene_taxon.region.abundance.txt \
	P22=P22/P22.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P42=P42/P42.gene_taxon.region.abundance.txt \
	> prediction.txt
