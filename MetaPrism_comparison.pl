#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use Statistics::R;
use Bio::DB::Taxonomy;

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

GetOptions(
	'h' => \(my $help = ''),
	'R=s' => \(my $taxonRank = 'genus'),
	'F=s' => \(my $featureType = 'gene_taxon'),
	't=s' => \(my $test = 'kruskal'),
	'o=f' => \(my $offset = 1),
);
my @featureTypeList = ('gene_taxon', 'gene', 'gene_average', 'gene_average_fraction', 'taxon', 'taxon_average', 'taxon_average_fraction');
my $featureTypes = join(', ', map {"\"$_\""} @featureTypeList);
my @testList = ('kruskal', 'anova', 'poisson', 'quasipoisson', 'metagenomeSeq');
my $tests = join(', ', map {"\"$_\""} @testList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_comparison.pl [options] sample.group.txt [sample=]abundance.txt [...] > MetaPrism_comparison.txt

Options: -h       display this help message
         -R STR   taxon rank [$taxonRank]
         -F STR   feature type, $featureTypes [$featureType]
         -t STR   statistical test for comparing sample groups, $tests [$test]
         -o FLOAT offset [$offset]

EOF
}
die "ERROR in $0: The feature type is not available.\n" if(scalar(grep {$featureType eq $_} @featureTypeList) == 0);
die "ERROR in $0: The test is not available.\n" if(scalar(grep {$test eq $_} @testList) == 0);

my ($groupFile, @sampleFileList) = @ARGV;

if(not -r "$dataPath/nodes.dmp" or not -r "$dataPath/names.dmp") {
	my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
	my $file = "$dataPath/taxdump.tar.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file);
	die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
	system("cd $dataPath; tar -zxf taxdump.tar.gz nodes.dmp");
	system("cd $dataPath; tar -zxf taxdump.tar.gz names.dmp");
	system("rm -f $dataPath/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
}
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataPath, -nodesfile => "$dataPath/nodes.dmp", -namesfile => "$dataPath/names.dmp");

my %sampleNumberListHash = ();
my %numberGroupHash = ();
my $number = 0;
{
	open(my $reader, $groupFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($sample, $group) = split(/\t/, $line, -1);
		$number += 1;
		push(@{$sampleNumberListHash{$sample}}, $number);
		$numberGroupHash{$number} = $group;
	}
	close($reader);
}

my %numberFeatureAbundanceHash = ();
my %featureHash = ();
my %geneDefinitionHash = ();
@sampleFileList = map {[$_->[0], $_->[-1]]} map {[split(/=/, $_, 2)]} @sampleFileList;
foreach(@sampleFileList) {
	my ($sample, $file) = @$_;
	if(-s $file && defined(my $numberList = $sampleNumberListHash{$sample})) {
		open(my $reader, $file);
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			my $taxonName = getTaxonName($tokenHash{'taxonId'});
			my $feature = '';
			if($featureType eq 'gene_taxon') {
				$featureHash{$feature = join("\t", $tokenHash{'gene'}, $taxonName)} = 1;
			} elsif($featureType eq 'gene' || $featureType eq 'gene_average' || $featureType eq 'gene_average_fraction') {
				$featureHash{$feature = $tokenHash{'gene'}}->{$taxonName} = 1;
			} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average' || $featureType eq 'taxon_average_fraction') {
				$featureHash{$feature = $taxonName}->{$tokenHash{'gene'}} = 1;
			}
			$numberFeatureAbundanceHash{$_}->{$feature} += $tokenHash{'abundance'} foreach(@$numberList);
			$geneDefinitionHash{$tokenHash{'gene'}} = $tokenHash{'definition'} if(defined($tokenHash{'definition'}) && $tokenHash{'definition'} ne '');
		}
		close($reader);
	}
}
if($featureType eq 'gene_average' || $featureType eq 'taxon_average' || $featureType eq 'gene_average_fraction' || $featureType eq 'taxon_average_fraction') {
	foreach my $feature (keys %featureHash) {
		my $count = scalar(keys %{$featureHash{$feature}});
		$numberFeatureAbundanceHash{$_}->{$feature} /= $count foreach(grep {defined($numberFeatureAbundanceHash{$_}->{$feature})} 1 .. $number);
	}
}
if($featureType eq 'gene_average_fraction' || $featureType eq 'taxon_average_fraction') {
	foreach my $number (1 .. $number) {
		my $sum = 0;
		$sum += $numberFeatureAbundanceHash{$number}->{$_} foreach(keys %{$numberFeatureAbundanceHash{$number}});
		$numberFeatureAbundanceHash{$number}->{$_} /= $sum foreach(keys %{$numberFeatureAbundanceHash{$number}});
	}
}

if(my @featureList = sort keys %featureHash) {
	my $R = Statistics::R->new();
	$R->run('x <- data.frame()');
	$R->run('y <- c()');
	foreach my $number (1 .. $number) {
		foreach my $index (0 .. $#featureList) {
			if(defined(my $abundance = $numberFeatureAbundanceHash{$number}->{$featureList[$index]})) {
				$R->set(sprintf('x[%d, %d]', $number, $index + 1), $abundance);
			} else {
				$R->set(sprintf('x[%d, %d]', $number, $index + 1), 0);
			}
		}
		$R->set("y[$number]", $numberGroupHash{$number});
	}
	foreach my $index (0 .. $#featureList) {
		$R->set(sprintf('colnames(x)[%d]', $index + 1), $featureList[$index]);
	}
	$R->run('x <- data.matrix(x)');
	$R->run('y <- factor(y)');

	$R->run(sprintf('table <- t(x + %f)', $offset));
	$R->run('condition <- y');
	$R->run('group <- apply(table, 1, function(x) {m <- tapply(x, condition, mean); n <- names(m)[max(m) == m]; paste(n, sep = ",", collapse = ",")})');
	$R->run('log2foldchange <- apply(table, 1, function(x) {m <- tapply(x, condition, mean); n <- names(m)[max(m) == m]; log2(mean(x[condition %in% n]) / mean(x[!(condition %in% n)]))})');
	if($test eq 'kruskal') {
		$R->run('p.value <- apply(table, 1, function(x) {kruskal.test(x, condition)$p.value})');
	}
	if($test eq 'anova') {
		$R->run('p.value <- apply(log2(table + 1), 1, function(x) {summary(aov(x ~ condition))[[1]]["condition", "Pr(>F)"]})');
	}
	if($test eq 'poisson') {
		$R->run('offset.log.colSums <- offset(log(colSums(table)))');
		$R->run('p.value <- apply(table, 1, function(x) {anova(glm(round(x) ~ 1 + offset.log.colSums + condition, family = "poisson"), test = "Chisq")["condition", "Pr(>Chi)"]})');
	}
	if($test eq 'quasipoisson') {
		$R->run('offset.log.colSums <- offset(log(colSums(table)))');
		$R->run('p.value <- apply(table, 1, function(x) {anova(glm(round(x) ~ 1 + offset.log.colSums + condition, family = "quasipoisson"), test = "Chisq")["condition", "Pr(>Chi)"]})');
	}
	if($test eq 'metagenomeSeq') {
		$R->run('library(metagenomeSeq)');
		$R->run('phenoData <- AnnotatedDataFrame(data.frame(condition = condition, row.names = colnames(table)))');
		$R->run('featureData <- AnnotatedDataFrame(data.frame(gene = rownames(table), row.names = rownames(table)))');
		$R->run('MRexperiment <- newMRexperiment(table, phenoData = phenoData, featureData = featureData)');
		$R->run('MRexperiment.p <- cumNormStat(MRexperiment, pFlag = TRUE)');
		$R->run('MRexperiment <- cumNorm(MRexperiment, p = MRexperiment.p)');
		$R->run('normFactor <- normFactors(MRexperiment)');
		$R->run('normFactor <- log2(normFactor / median(normFactor) + 1)');
		$R->run('mod <- model.matrix(~ condition + normFactor)');
		$R->run('fit <- fitZig(obj = MRexperiment, mod = mod)');
		$R->run('MRcoefs <- MRcoefs(fit, number = nrow(assayData(MRexperiment)$counts))');
		$R->run('table <- MRcounts(MRexperiment, norm = TRUE, log = TRUE)[rownames(table), colnames(table)]');
		$R->run('p.value <- MRcoefs[rownames(table), "pvalues"]');
		$R->run('p.value[is.na(p.value)] <- 1');
	}
	$R->run('p.adjust <- p.adjust(p.value, method = "fdr")');
	die "ERROR in $0: R.\n" if(grep {$R->get(sprintf('length(%s)', $_)) != scalar(@featureList)} ('group', 'log2foldchange', 'p.value', 'p.adjust'));

	if($featureType eq 'gene_taxon') {
		if(%geneDefinitionHash) {
			my @columnList = ('gene', 'taxon', 'definition', 'group', 'log2foldchange', 'p.value', 'p.adjust');
			print join("\t", @columnList), "\n";
			foreach my $index (0 .. $#featureList) {
				my %tokenHash = ();
				@tokenHash{'gene', 'taxon'} = split(/\t/, $featureList[$index]);
				$tokenHash{'definition'} = $geneDefinitionHash{$tokenHash{'gene'}};
				@tokenHash{'group', 'log2foldchange', 'p.value', 'p.adjust'} = map {$R->get(sprintf('%s[%d]', $_, $index + 1))} ('group', 'log2foldchange', 'p.value', 'p.adjust');
				print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		} else {
			my @columnList = ('gene', 'taxon', 'group', 'log2foldchange', 'p.value', 'p.adjust');
			print join("\t", @columnList), "\n";
			foreach my $index (0 .. $#featureList) {
				my %tokenHash = ();
				@tokenHash{'gene', 'taxon'} = split(/\t/, $featureList[$index]);
				@tokenHash{'group', 'log2foldchange', 'p.value', 'p.adjust'} = map {$R->get(sprintf('%s[%d]', $_, $index + 1))} ('group', 'log2foldchange', 'p.value', 'p.adjust');
				print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		}
	} elsif($featureType eq 'gene' || $featureType eq 'gene_average' || $featureType eq 'gene_average_fraction') {
		if(%geneDefinitionHash) {
			my @columnList = ('gene', 'definition', 'group', 'log2foldchange', 'p.value', 'p.adjust');
			print join("\t", @columnList), "\n";
			foreach my $index (0 .. $#featureList) {
				my %tokenHash = ();
				$tokenHash{'gene'} = $featureList[$index];
				$tokenHash{'definition'} = $geneDefinitionHash{$tokenHash{'gene'}};
				@tokenHash{'group', 'log2foldchange', 'p.value', 'p.adjust'} = map {$R->get(sprintf('%s[%d]', $_, $index + 1))} ('group', 'log2foldchange', 'p.value', 'p.adjust');
				print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		} else {
			my @columnList = ('gene', 'group', 'log2foldchange', 'p.value', 'p.adjust');
			print join("\t", @columnList), "\n";
			foreach my $index (0 .. $#featureList) {
				my %tokenHash = ();
				$tokenHash{'gene'} = $featureList[$index];
				@tokenHash{'group', 'log2foldchange', 'p.value', 'p.adjust'} = map {$R->get(sprintf('%s[%d]', $_, $index + 1))} ('group', 'log2foldchange', 'p.value', 'p.adjust');
				print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		}
	} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average' || $featureType eq 'taxon_average_fraction') {
		my @columnList = ('taxon', 'group', 'log2foldchange', 'p.value', 'p.adjust');
		print join("\t", @columnList), "\n";
		foreach my $index (0 .. $#featureList) {
			my %tokenHash = ();
			$tokenHash{'taxon'} = $featureList[$index];
			@tokenHash{'group', 'log2foldchange', 'p.value', 'p.adjust'} = map {$R->get(sprintf('%s[%d]', $_, $index + 1))} ('group', 'log2foldchange', 'p.value', 'p.adjust');
			print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
		}
	}
	$R->stop();
}

sub getTaxonName {
	my ($taxonId) = @_;
	if(defined($taxonId) && defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
		do {
			return $taxon->scientific_name if($taxon->rank eq $taxonRank);
		} while(defined($taxon = $taxon->ancestor));
	}
	return 'undefined';
}
