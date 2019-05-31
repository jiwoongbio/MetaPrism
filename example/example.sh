#!/bin/env bash

# Prepare database files
perl ../MetaPrism_gene_prepare.pl

# Gene-taxon annotation and abundance quantification
#(cd P1; ./P1.sh)
#(cd P2; ./P2.sh)
#(cd P3; ./P3.sh)
#(cd P4; ./P4.sh)
#(cd P5; ./P5.sh)
#(cd P6; ./P6.sh)
#(cd P7; ./P7.sh)
#(cd P8; ./P8.sh)
#(cd P9; ./P9.sh)
#(cd P10; ./P10.sh)
#(cd P12; ./P12.sh)
#(cd P19; ./P19.sh)
#(cd P31; ./P31.sh)
#(cd P55; ./P55.sh)
#(cd P63; ./P63.sh)
#(cd P78; ./P78.sh)

# Compare sample groups and identify differentially-abundant genes
perl ../MetaPrism_comparison.pl -F gene sample.group.txt \
	P1=P1/P1.gene_taxon.region.abundance.txt \
	P2=P2/P2.gene_taxon.region.abundance.txt \
	P3=P3/P3.gene_taxon.region.abundance.txt \
	P4=P4/P4.gene_taxon.region.abundance.txt \
	P5=P5/P5.gene_taxon.region.abundance.txt \
	P6=P6/P6.gene_taxon.region.abundance.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P9=P9/P9.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P12=P12/P12.gene_taxon.region.abundance.txt \
	P19=P19/P19.gene_taxon.region.abundance.txt \
	P31=P31/P31.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P78=P78/P78.gene_taxon.region.abundance.txt \
	> gene.comparison.txt

# Prediction model accuracy
perl ../MetaPrism_prediction.pl sample.group.txt \
	P1=P1/P1.gene_taxon.region.abundance.txt \
	P2=P2/P2.gene_taxon.region.abundance.txt \
	P3=P3/P3.gene_taxon.region.abundance.txt \
	P4=P4/P4.gene_taxon.region.abundance.txt \
	P5=P5/P5.gene_taxon.region.abundance.txt \
	P6=P6/P6.gene_taxon.region.abundance.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P9=P9/P9.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P12=P12/P12.gene_taxon.region.abundance.txt \
	P19=P19/P19.gene_taxon.region.abundance.txt \
	P31=P31/P31.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P78=P78/P78.gene_taxon.region.abundance.txt \
	> prediction.txt

# Generate table
perl ../MetaPrism_table.pl -F taxon_average -s \
	P1=P1/P1.gene_taxon.region.abundance.txt \
	P2=P2/P2.gene_taxon.region.abundance.txt \
	P3=P3/P3.gene_taxon.region.abundance.txt \
	P4=P4/P4.gene_taxon.region.abundance.txt \
	P5=P5/P5.gene_taxon.region.abundance.txt \
	P6=P6/P6.gene_taxon.region.abundance.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P9=P9/P9.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P12=P12/P12.gene_taxon.region.abundance.txt \
	P19=P19/P19.gene_taxon.region.abundance.txt \
	P31=P31/P31.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P78=P78/P78.gene_taxon.region.abundance.txt \
	> taxon.table.txt

awk -F'\t' '(NR == 1 || ($4 > 1 && $5 < 0.01))' gene.comparison.txt > gene.comparison.filtered.txt

# Generate heatmap
perl ../MetaPrism_heatmap.pl -F gene -s -g gene.comparison.filtered.txt -r both \
	P1=P1/P1.gene_taxon.region.abundance.txt \
	P2=P2/P2.gene_taxon.region.abundance.txt \
	P3=P3/P3.gene_taxon.region.abundance.txt \
	P4=P4/P4.gene_taxon.region.abundance.txt \
	P5=P5/P5.gene_taxon.region.abundance.txt \
	P6=P6/P6.gene_taxon.region.abundance.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P9=P9/P9.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P12=P12/P12.gene_taxon.region.abundance.txt \
	P19=P19/P19.gene_taxon.region.abundance.txt \
	P31=P31/P31.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P78=P78/P78.gene_taxon.region.abundance.txt \
	> gene.heatmap.html
