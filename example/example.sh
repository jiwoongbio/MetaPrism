#!/bin/env bash

# Prepare database files
perl ../MetaPrism_gene_prepare.pl

# Gene-taxon annotation and abundance quantification
(cd P7; ./P7.sh)
(cd P8; ./P8.sh)
(cd P10; ./P10.sh)
(cd P14; ./P14.sh)
(cd P23; ./P23.sh)
(cd P25; ./P25.sh)
(cd P34; ./P34.sh)
(cd P35; ./P35.sh)
(cd P52; ./P52.sh)
(cd P55; ./P55.sh)
(cd P59; ./P59.sh)
(cd P61; ./P61.sh)
(cd P63; ./P63.sh)
(cd P67; ./P67.sh)
(cd P68; ./P68.sh)
(cd P69; ./P69.sh)
(cd P24; ./P24.sh)
(cd P28; ./P28.sh)
(cd P30; ./P30.sh)
(cd P32; ./P32.sh)
(cd P46; ./P46.sh)
(cd P53; ./P53.sh)
(cd P54; ./P54.sh)
(cd P66; ./P66.sh)

# Compare sample groups and identify differentially-abundant genes
perl ../MetaPrism_comparison.pl -F gene sample.group.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P35=P35/P35.gene_taxon.region.abundance.txt \
	P52=P52/P52.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P59=P59/P59.gene_taxon.region.abundance.txt \
	P61=P61/P61.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P67=P67/P67.gene_taxon.region.abundance.txt \
	P68=P68/P68.gene_taxon.region.abundance.txt \
	P69=P69/P69.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P28=P28/P28.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P46=P46/P46.gene_taxon.region.abundance.txt \
	P53=P53/P53.gene_taxon.region.abundance.txt \
	P54=P54/P54.gene_taxon.region.abundance.txt \
	P66=P66/P66.gene_taxon.region.abundance.txt \
	> gene.comparison.txt

awk -F'\t' '(NR == 1 || ($4 >= 1 && $5 <= 0.01))' gene.comparison.txt > gene.comparison.filtered.txt

# Generate heatmap
perl ../MetaPrism_heatmap.pl -F gene -s -g gene.comparison.filtered.txt -r both \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P35=P35/P35.gene_taxon.region.abundance.txt \
	P52=P52/P52.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P59=P59/P59.gene_taxon.region.abundance.txt \
	P61=P61/P61.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P67=P67/P67.gene_taxon.region.abundance.txt \
	P68=P68/P68.gene_taxon.region.abundance.txt \
	P69=P69/P69.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P28=P28/P28.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P46=P46/P46.gene_taxon.region.abundance.txt \
	P53=P53/P53.gene_taxon.region.abundance.txt \
	P54=P54/P54.gene_taxon.region.abundance.txt \
	P66=P66/P66.gene_taxon.region.abundance.txt \
	> gene.heatmap.html

# Generate table
perl ../MetaPrism_table.pl -F taxon_average -s \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P35=P35/P35.gene_taxon.region.abundance.txt \
	P52=P52/P52.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P59=P59/P59.gene_taxon.region.abundance.txt \
	P61=P61/P61.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P67=P67/P67.gene_taxon.region.abundance.txt \
	P68=P68/P68.gene_taxon.region.abundance.txt \
	P69=P69/P69.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P28=P28/P28.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P46=P46/P46.gene_taxon.region.abundance.txt \
	P53=P53/P53.gene_taxon.region.abundance.txt \
	P54=P54/P54.gene_taxon.region.abundance.txt \
	P66=P66/P66.gene_taxon.region.abundance.txt \
	> taxon.table.txt

# Prediction model accuracy
perl ../MetaPrism_prediction.pl -t xgbTree -f prediction.feature.txt sample.group.txt \
	P7=P7/P7.gene_taxon.region.abundance.txt \
	P8=P8/P8.gene_taxon.region.abundance.txt \
	P10=P10/P10.gene_taxon.region.abundance.txt \
	P14=P14/P14.gene_taxon.region.abundance.txt \
	P23=P23/P23.gene_taxon.region.abundance.txt \
	P25=P25/P25.gene_taxon.region.abundance.txt \
	P34=P34/P34.gene_taxon.region.abundance.txt \
	P35=P35/P35.gene_taxon.region.abundance.txt \
	P52=P52/P52.gene_taxon.region.abundance.txt \
	P55=P55/P55.gene_taxon.region.abundance.txt \
	P59=P59/P59.gene_taxon.region.abundance.txt \
	P61=P61/P61.gene_taxon.region.abundance.txt \
	P63=P63/P63.gene_taxon.region.abundance.txt \
	P67=P67/P67.gene_taxon.region.abundance.txt \
	P68=P68/P68.gene_taxon.region.abundance.txt \
	P69=P69/P69.gene_taxon.region.abundance.txt \
	P24=P24/P24.gene_taxon.region.abundance.txt \
	P28=P28/P28.gene_taxon.region.abundance.txt \
	P30=P30/P30.gene_taxon.region.abundance.txt \
	P32=P32/P32.gene_taxon.region.abundance.txt \
	P46=P46/P46.gene_taxon.region.abundance.txt \
	P53=P53/P53.gene_taxon.region.abundance.txt \
	P54=P54/P54.gene_taxon.region.abundance.txt \
	P66=P66/P66.gene_taxon.region.abundance.txt \
	> prediction.txt
