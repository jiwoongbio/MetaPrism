# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use Statistics::R;
use Bio::DB::Taxonomy;

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/MetaPrism_data";
system("mkdir -p $dataPath");

GetOptions(
	'h' => \(my $help = ''),
	'A=s' => \(my $abundanceColumn = 'meanDepth/genome'),
	'R=s' => \(my $taxonRank = 'genus'),
	'F=s' => \(my $featureType = 'gene_taxon'),
	't=s' => \(my $trainMethod = 'rf'),
	'c=s' => \(my $trainControlMethod = 'cv'),
	'm=s' => \(my $modelFile = ''),
	'f=s' => \(my $featureImportanceFile = ''),
	's=i' => \(my $seed = 1),
);
my @featureTypeList = ('gene_taxon', 'gene', 'gene_average', 'gene_average_fraction', 'taxon', 'taxon_average', 'taxon_average_fraction');
my $featureTypes = join(', ', map {"\"$_\""} @featureTypeList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_prediction.pl [options] sample.group.txt [sample=]abundance.txt [...]

Options: -h       display this help message
         -A STR   abundance column [$abundanceColumn]
         -R STR   taxon rank [$taxonRank]
         -F STR   feature type, $featureTypes [$featureType]
         -t STR   train method [$trainMethod]
         -c STR   train control method [$trainControlMethod]
         -m FILE  model file
         -f FILE  important feature file
         -s INT   seed [$seed]

EOF
}
die "ERROR in $0: The feature type is not available.\n" if(scalar(grep {$featureType eq $_} @featureTypeList) == 0);

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
			my $taxonName = getTaxonName($tokenHash{'taxon'});
			my $feature = '';
			if($featureType eq 'gene_taxon') {
				$featureHash{$feature = join("\t", $tokenHash{'gene'}, $taxonName)} = 1;
			} elsif($featureType eq 'gene' || $featureType eq 'gene_average' || $featureType eq 'gene_average_fraction') {
				$featureHash{$feature = $tokenHash{'gene'}}->{$taxonName} = 1;
			} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average' || $featureType eq 'taxon_average_fraction') {
				$featureHash{$feature = $taxonName}->{$tokenHash{'gene'}} = 1;
			}
			$numberFeatureAbundanceHash{$_}->{$feature} += $tokenHash{$abundanceColumn} foreach(@$numberList);
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

	$R->run('library(caret)');
	$R->run(sprintf('set.seed(%d)', $seed));
	$R->run(sprintf('model <- train(x, y, %s)', join(', ',
		sprintf('method = "%s"', $trainMethod),
		sprintf('trControl = trainControl(method = "%s")', $trainControlMethod),
	)));
	$R->run(sprintf('save(model, file = "%s")', $modelFile)) if($modelFile ne '');
	$R->run(sprintf('write.table(varImp(model)$importance, file = "%s", quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)', $featureImportanceFile)) if($featureImportanceFile ne '');
	{
		my $ncol = $R->get('ncol(model$results)');
		my $nrow = $R->get('nrow(model$results)');
		print join("\t", map {$R->get(sprintf('colnames(model$results)[%d]', $_))} 1 .. $ncol), "\n";
		foreach my $row (1 .. $nrow) {
			print join("\t", map {$R->get(sprintf('model$results[%d, %d]', $row, $_))} 1 .. $ncol), "\n";
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
