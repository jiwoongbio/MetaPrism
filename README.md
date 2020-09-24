# MetaPrism
Taxonomical and Functional Profiling for Shotgun Metagenomic Sequencing


## Requirements

1. Perl - https://www.perl.org
2. R - http://www.r-project.org
3. Perl module Statistics::R - https://metacpan.org/pod/Statistics::R
4. R library caret - https://cran.r-project.org/web/packages/caret/index.html
5. R library randomForest
6. DIAMOND - https://github.com/bbuchfink/diamond or USEARCH - https://www.drive5.com/usearch/
7. BWA - http://bio-bwa.sourceforge.net
8. Samtools - http://www.htslib.org
9. Centrifuge - https://ccb.jhu.edu/software/centrifuge/
10. Linux commands: sort, wget - https://www.gnu.org/software/wget/


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git. It will take a few seconds.
```
git clone https://github.com/jiwoongbio/MetaPrism.git
```


## Usages

* MetaPrism_gene_prepare.pl
```
Usage:   perl MetaPrism_gene_preapare.pl [options]

Options: -h       display this help message
         -r       redownload data
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -k       download prebuilt KEGG files
         -a       download ARDB database
         -b       download beta-lactamase database
```

* MetaPrism_gene.pl
```
Usage:   perl MetaPrism_gene.pl [options] output.prefix genome.fasta [input.fastq|input.R1.fastq,input.R2.fastq [...]]

Options: -h       display this help message
         -A STR   prepared genome prefix
         -B       input indexed sorted BAM file instead of FASTQ file
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments [10]
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -a FLOAT search acceleration for ublast [0.5]
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [GTG,ATG,CTG,TTG,ATA,ATC,ATT]
         -T STR   comma-separated termination codons [TAG,TAA,TGA]
         -l INT   minimum translation length [10]
         -c FLOAT minimum coverage [0.8]
         -q INT   minimum mapping quality [0]
         -s STR   strand specificity, "f" or "r"
         -P STR   contig prefix used for abundance estimation
```

* MetaPrism_taxon_centrifuge.pl
```
Usage:   perl MetaPrism_taxon_centrifuge.pl [options] MetaPrism_gene.region.txt genome.fasta centrifuge.index > MetaPrism_gene.region.taxon.txt

Options: -h       display this help message
         -p INT   number of threads [1]
```

* MetaPrism_comparison.pl
```
Usage:   perl MetaPrism_comparison.pl [options] sample.group.txt [sample=]abundance.txt [...] > MetaPrism_comparison.txt

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -t STR   statistical test for comparing sample groups, "kruskal", "anova", "poisson", "quasipoisson", "metagenomeSeq" [kruskal]
         -o FLOAT offset [1]
```

* MetaPrism_prediction.pl
```
Usage:   perl MetaPrism_prediction.pl [options] sample.group.txt [sample=]abundance.txt [...]

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -t STR   train method [rf]
         -c STR   train control method [LOOCV]
         -m FILE  model file
         -f FILE  important feature file
         -s INT   seed [1]
```

* MetaPrism_table.pl
```
Usage:   perl MetaPrism_table.pl [options] [sample=]abundance.txt [...] > table.txt

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -s       scale
```

* MetaPrism_heatmap.pl
```
Usage:   perl MetaPrism_heatmap.pl [options] [sample=]abundance.txt [...] > heatmap.html

Options: -h       display this help message
         -A STR   abundance column [meanDepth/genome]
         -R STR   taxon rank [genus]
         -F STR   feature type, "gene_taxon", "gene", "gene_average", "taxon", "taxon_average" [gene_taxon]
         -s       scale
         -g FILE  feature file
         -t INT   taxon abbreviation length [4]
         -f INT   HTML font size [15]
         -w INT   HTML table cell width [60]
```
