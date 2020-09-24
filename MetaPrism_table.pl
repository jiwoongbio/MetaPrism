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
	's' => \(my $scale = ''),
);
my @featureTypeList = ('gene_taxon', 'gene', 'gene_average', 'taxon', 'taxon_average');
my $featureTypes = join(', ', map {"\"$_\""} @featureTypeList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_table.pl [options] [sample=]abundance.txt [...] > table.txt

Options: -h       display this help message
         -A STR   abundance column [$abundanceColumn]
         -R STR   taxon rank [$taxonRank]
         -F STR   feature type, $featureTypes [$featureType]
         -s       scale

EOF
}
die "ERROR in $0: The feature type is not available.\n" if(scalar(grep {$featureType eq $_} @featureTypeList) == 0);

my (@sampleFileList) = @ARGV;

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

my @sampleList = ();
my %sampleFeatureAbundanceHash = ();
my %sampleFeatureCountHash = ();
my %geneDefinitionHash = ();
@sampleFileList = map {[$_->[0], $_->[-1]]} map {[split(/=/, $_, 2)]} @sampleFileList;
foreach(@sampleFileList) {
	my ($sample, $file) = @$_;
	if(-s $file) {
		push(@sampleList, $sample);
		open(my $reader, $file);
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			my $feature = '';
			if($featureType eq 'gene_taxon') {
				if((my $taxonName = getTaxonName($tokenHash{'taxon'})) ne '') {
					$feature = join("\t", $tokenHash{'gene'}, $taxonName);
				}
			} elsif($featureType eq 'gene' || $featureType eq 'gene_average') {
				$feature = $tokenHash{'gene'};
			} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average') {
				if((my $taxonName = getTaxonName($tokenHash{'taxon'})) ne '') {
					$feature = $taxonName;
				}
			}
			if($feature ne '') {
				$sampleFeatureAbundanceHash{$sample}->{$feature} += $tokenHash{$abundanceColumn};
				$sampleFeatureCountHash{$sample}->{$feature} += 1;
			}
			$geneDefinitionHash{$tokenHash{'gene'}} = $tokenHash{'definition'} if(defined($tokenHash{'definition'}) && $tokenHash{'definition'} ne '');
		}
		close($reader);
	}
}
if($featureType eq 'gene_average' || $featureType eq 'taxon_average') {
	foreach my $sample (@sampleList) {
		$sampleFeatureAbundanceHash{$sample}->{$_} /= $sampleFeatureCountHash{$sample}->{$_} foreach(keys %{$sampleFeatureAbundanceHash{$sample}});
	}
}
my %featureHash = ();
foreach my $sample (@sampleList) {
	$featureHash{$_} = 1 foreach(keys %{$sampleFeatureAbundanceHash{$sample}});
}

if(my @featureList = sort keys %featureHash) {
	if($featureType eq 'gene_taxon') {
		if(%geneDefinitionHash) {
			my @featureColumnList = ('gene', 'taxon', 'definition');
			print join("\t", @featureColumnList, @sampleList), "\n";
			foreach my $feature (@featureList) {
				my %tokenHash = ();
				@tokenHash{'gene', 'taxon'} = split(/\t/, $feature);
				$tokenHash{'definition'} = $geneDefinitionHash{$tokenHash{'gene'}};
				my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
				@abundanceList = scale(@abundanceList) if($scale);
				print join("\t", (map {defined($_) ? $_ : ''} @tokenHash{@featureColumnList}), @abundanceList), "\n";
			}
		} else {
			my @featureColumnList = ('gene', 'taxon');
			print join("\t", @featureColumnList, @sampleList), "\n";
			foreach my $feature (@featureList) {
				my %tokenHash = ();
				@tokenHash{'gene', 'taxon'} = split(/\t/, $feature);
				my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
				@abundanceList = scale(@abundanceList) if($scale);
				print join("\t", (map {defined($_) ? $_ : ''} @tokenHash{@featureColumnList}), @abundanceList), "\n";
			}
		}
	} elsif($featureType eq 'gene' || $featureType eq 'gene_average') {
		if(%geneDefinitionHash) {
			my @featureColumnList = ('gene', 'definition');
			print join("\t", @featureColumnList, @sampleList), "\n";
			foreach my $feature (@featureList) {
				my %tokenHash = ();
				$tokenHash{'gene'} = $feature;
				$tokenHash{'definition'} = $geneDefinitionHash{$tokenHash{'gene'}};
				my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
				@abundanceList = scale(@abundanceList) if($scale);
				print join("\t", (map {defined($_) ? $_ : ''} @tokenHash{@featureColumnList}), @abundanceList), "\n";
			}
		} else {
			my @featureColumnList = ('gene');
			print join("\t", @featureColumnList, @sampleList), "\n";
			foreach my $feature (@featureList) {
				my %tokenHash = ();
				$tokenHash{'gene'} = $feature;
				my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
				@abundanceList = scale(@abundanceList) if($scale);
				print join("\t", (map {defined($_) ? $_ : ''} @tokenHash{@featureColumnList}), @abundanceList), "\n";
			}
		}
	} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average') {
		my @featureColumnList = ('taxon');
		print join("\t", @featureColumnList, @sampleList), "\n";
		foreach my $feature (@featureList) {
			my %tokenHash = ();
			$tokenHash{'taxon'} = $feature;
			my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
			@abundanceList = scale(@abundanceList) if($scale);
			print join("\t", (map {defined($_) ? $_ : ''} @tokenHash{@featureColumnList}), @abundanceList), "\n";
		}
	}
}

sub scale {
	my @valueList = @_;
	my $mean = 0;
	$mean += $_ foreach(@valueList);
	$mean /= scalar(@valueList);
	my $sd = 0;
	$sd += ($_ - $mean) ** 2 foreach(@valueList);
	$sd /= scalar(@valueList) - 1;
	$sd = $sd ** 0.5;
	return map {($_ - $mean) / $sd} @valueList;
}

sub getTaxonName {
	my ($taxonId) = @_;
	if(defined(my $taxon = $db->get_taxon(-taxonid => $taxonId))) {
		do {
			return $taxon->scientific_name if($taxon->rank eq $taxonRank);
		} while(defined($taxon = $taxon->ancestor));
	}
	return '';
}
