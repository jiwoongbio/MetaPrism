Table of Contents
- [MetaPrism](#metaprism)
- [Requirements](#requirements)
- [Install](#install)
- [Tutorial](#tutorial)
- [Usages](#usages)
- [Citation](#citation)

# MetaPrism
MetaPrism: A Toolkit for Joint Analysis of Meta-genomic Sequencing Data 

MetaPrism provides joint profile (infer both taxonomical and functional profile) for shotgun metagenomic sequencing data. It also offer tools to 1) classify sequence reads and estimate the abundances for taxa-specific genes; 2) tabularize and visualize taxa-specific gene abundances; 3) build asso-ciation and prediction models for comparative analysis. 


## Requirements

1. Perl - https://www.perl.org
2. R - http://www.r-project.org
3. Perl module Statistics::R - https://metacpan.org/pod/Statistics::R
4. R library caret - https://cran.r-project.org/web/packages/caret/index.html
5. R library randomForest - https://cran.r-project.org/web/packages/randomForest/index.html
6. DIAMOND - https://github.com/bbuchfink/diamond (recommended) or USEARCH - https://www.drive5.com/usearch/
7. BWA - http://bio-bwa.sourceforge.net
8. Samtools - http://www.htslib.org
9. Centrifuge - https://ccb.jhu.edu/software/centrifuge/
10. Linux commands: sort, wget - https://www.gnu.org/software/wget/
11. MEGAHIT - https://github.com/voutcn/megahit (optional) 


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git. It will take a few seconds.
```
git clone https://github.com/jiwoongbio/MetaPrism.git
```


## Tutorial

We present a short tutorial to help users quickly get started on their own analysis. The datasets are based on [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397906](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397906) and the full command list is  available at [example/example.sh](example/example.sh).

1. Prepare database files
This step will download necessary databases for MetaPrism in your current directory.
```
perl MetaPrism_gene_prepare.pl
```

2. De novo metagenome assembly (per sample)
For each sample, users need to perform de novo assembly. Suppose the input paired end sequence files are `sample1.1.fastq.gz` and `sample1.2.fastq.gz`. The commnad line is:
```
megahit -1 sample1.1.fastq.gz -2 sample1.2.fastq.gz -o sample1.megahit
```

3. Gene annotation and abundance quantification (per sample)
For each sample, MetaPrism first quantify all gene abundances. The file `sample1.megahit/final.contigs.fa` is the output from the previous step.
```
perl MetaPrism_gene.pl sample1.gene sample1.megahit/final.contigs.fa sample1.1.fastq.gz,sample1.2.fastq.gz
```

4. Taxon annotation (per sample)
MetaPrism will next infer the taxonomy for each contig. The file `sample1.gene.region.abundance.txt` is the output from the previous step. The result will be outputted to the console. We redirected it to the result file `sample1.gene_taxon.region.abundance.txt`.
```
perl MetaPrism_taxon_centrifuge.pl sample1.gene.region.abundance.txt sample1.megahit/final.contigs.fa centrifuge/data/p_compressed > sample1.gene_taxon.region.abundance.txt
```

5. Compare sample groups and identify differentially-abundant genes
Suppose that you repeat step 2 to step 4, and you get a list of joint features (`sample1.gene_taxon.region.abundance.txt`, `sample2.gene_taxon.region.abundance.txt`, ..., `sample6.gene_taxon.region.abundance.txt`). MetaPrism can perform comparative analysis using the following command:
```
perl MetaPrism_comparison.pl -F gene sample.group.txt \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> gene.comparison.txt

awk -F'\t' '(NR == 1 || ($4 >= 1 && $5 <= 0.01))' gene.comparison.txt > gene.comparison.filtered.txt
```

* sample.group.txt is a text file containing lines of tab-delimited sample and group like following:

````
  sample1 group1
  sample2 group1
  sample3 group1
  sample4 group2
  sample5 group2
  sample6 group2
````

Here `sample1`, `sample2`, and `sample3` are from `group1`, and the rest are from `group2`.
For another example group file, see [example/sample.group.txt](example/sample.group.txt).

6. Generate a heatmap webpage
You can also generate a heatmap webpage using the `MetaPrism_heatmap.pl` command.
```
perl MetaPrism_heatmap.pl -F gene -s -g gene.comparison.filtered.txt -r both \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> gene.heatmap.html
```

7. Generate a tabular result file
Users may also want to have a tabular file for their own analysis. This command will produce such tabular text file: 
```
perl MetaPrism_table.pl -F taxon_average -s \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> taxon.table.txt
```

8. Build a prediction model
Users can build a prediction model using `MetaPrism_prediction.pl`. The option `-t xgbTree` uses the xgboost algorithm with leave-one-out cross validation. The result file `prediction.feature.txt` lists the feature importances, and the file `prediction.txt` lists the prediction accuracies.
```
perl MetaPrism_prediction.pl -t xgbTree -f prediction.feature.txt sample.group.txt \
	sample1=sample1.gene_taxon.region.abundance.txt \
	sample2=sample2.gene_taxon.region.abundance.txt \
	sample3=sample3.gene_taxon.region.abundance.txt \
	sample4=sample4.gene_taxon.region.abundance.txt \
	sample5=sample5.gene_taxon.region.abundance.txt \
	sample6=sample6.gene_taxon.region.abundance.txt \
	> prediction.txt
```


## Usages

In this section, we list the command line option for all available `MetaPrism` functions.
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


## Citation

MetaPrism: A Toolkit for Joint Taxa/Gene Analysis of Metagenomic Sequencing Data
https://www.biorxiv.org/content/10.1101/664748v1
