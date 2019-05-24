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
	'g=s' => \(my $featureFile = ''),
	't=i' => \(my $taxonAbbreviationLength = 4),
	'f=i' => \(my $fontSize = 15),
	'w=i' => \(my $tdWidth = 60),
);
my @featureTypeList = ('gene_taxon', 'gene', 'gene_average', 'taxon', 'taxon_average');
my $featureTypes = join(', ', map {"\"$_\""} @featureTypeList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_heatmap.pl [options] [sample=]abundance.txt [...] > heatmap.html

Options: -h       display this help message
         -A STR   abundance column [$abundanceColumn]
         -R STR   taxon rank [$taxonRank]
         -F STR   feature type, $featureTypes [$featureType]
         -s       scale
         -g FILE  feature file
         -t INT   taxon abbreviation length [$taxonAbbreviationLength]
         -f INT   HTML font size [$fontSize]
         -w INT   HTML table cell width [$tdWidth]

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
my %sampleGeneTaxonAbundanceHash = ();
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
				if((my $taxonName = getTaxonName($tokenHash{'taxon'})) ne '') {
					$sampleGeneTaxonAbundanceHash{$sample}->{$feature}->{$taxonName} += $tokenHash{$abundanceColumn};
				}
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
my @featureList = ();
if($featureFile ne '') {
	die "ERROR in $0: '$featureFile' is not readable.\n" unless(-r $featureFile);
	open(my $reader, $featureFile);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my $feature = $featureType eq 'gene_taxon' ? join("\t", @tokenList[0, 1]) : $tokenList[0];
		push(@featureList, $feature) if($featureHash{$feature});
	}
	close($reader);
} else {
	@featureList = sort keys %featureHash;
	if(scalar(@featureList) <= 65536) {
		my $R = Statistics::R->new();
		$R->run('x <- data.frame()');
		foreach my $featureIndex (0 .. $#featureList) {
			foreach my $sampleIndex (0 .. $#sampleList) {
				if(defined(my $abundance = $sampleFeatureAbundanceHash{$sampleList[$sampleIndex]}->{$featureList[$featureIndex]})) {
					$R->set(sprintf('x[%d, %d]', $sampleIndex + 1, $featureIndex + 1), $abundance);
				} else {
					$R->set(sprintf('x[%d, %d]', $sampleIndex + 1, $featureIndex + 1), 0);
				}
			}
		}
		foreach my $sampleIndex (0 .. $#sampleList) {
			$R->set(sprintf('rownames(x)[%d]', $sampleIndex + 1), $sampleList[$sampleIndex]);
		}
		foreach my $featureIndex (0 .. $#featureList) {
			$R->set(sprintf('colnames(x)[%d]', $featureIndex + 1), $featureList[$featureIndex]);
		}
		$R->run('hc <- hclust(as.dist(1 - cor(x, method = "spearman")))');
		@featureList = @{$R->get('hc$labels[hc$order]')};
		$R->stop();
	}
}

if(@featureList) {
	print <<EOF;
<!DOCTYPE html>
<html>
<head>
<style>
* {
	font-size: ${fontSize}px;
	line-height: ${fontSize}px;
}
table {
	border-collapse: collapse;
	border-spacing: 0px;
}
td {
	vertical-align: top;
	white-space: nowrap;
	padding: 5px;
}
table.heatmap td {
	border: 1px solid #000000;
}
table.heatmap td div.popup {
	display: none;
	position: absolute;
	background-color: #FFFFFF;
	padding: 5px;
	border: 1px solid #000000;
}
table.heatmap td:hover div.popup {
	display: block;
}
</style>
</head>
<body>
<form method="post">
EOF

	print "<table><tr>\n";
	print "<td><table class=\"heatmap\">\n";
	{
		my @tdList = ();
		if($featureType eq 'gene_taxon') {
			push(@tdList, map {"<td>$_</td>"} 'gene', 'taxon');
		} elsif($featureType eq 'gene' || $featureType eq 'gene_average') {
			push(@tdList, map {"<td>$_</td>"} 'gene');
		} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average') {
			push(@tdList, map {"<td>$_</td>"} 'taxon');
		}
		foreach my $sample (@sampleList) {
			push(@tdList, "<td style=\"text-align: center; width: ${tdWidth}px;\">$sample</td>");
		}
		print join('', "<tr>", @tdList, "</tr>"), "\n";
	}
	my $maximumScale = 0;
	if($scale) {
		foreach my $feature (@featureList) {
			my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
			foreach(map {abs($_)} scale(@abundanceList)) {
				$maximumScale = $_ if($_ > $maximumScale);
			}
		}
		$maximumScale = int($maximumScale) == $maximumScale ? $maximumScale : int($maximumScale + 1);
	}
	foreach my $feature (@featureList) {
		my @tdList = ();
		if($featureType eq 'gene_taxon') {
			my ($gene, $taxon) = split(/\t/, $feature);
			if(defined(my $definition = $geneDefinitionHash{$gene})) {
				$gene = "$gene<div class=\"popup\">$definition</div>";
			}
			push(@tdList, map {"<td>$_</td>"} $gene, $taxon);
		} elsif($featureType eq 'gene' || $featureType eq 'gene_average') {
			my $gene = $feature;
			if(defined(my $definition = $geneDefinitionHash{$gene})) {
				$gene = "$gene<div class=\"popup\">$definition</div>";
			}
			push(@tdList, map {"<td>$_</td>"} $gene);
		} elsif($featureType eq 'taxon' || $featureType eq 'taxon_average') {
			my $taxon = $feature;
			push(@tdList, map {"<td>$_</td>"} $taxon);
		}
		my @abundanceList = map {defined($_) ? $_ : 0} map {$_->{$feature}} @sampleFeatureAbundanceHash{@sampleList};
		my @colorList = ();
		if($scale) {
			@colorList = map {$_ > 0 ? color(1, ($maximumScale - $_) / $maximumScale, 0) : color(($maximumScale - abs($_)) / $maximumScale, 1, 0)} scale(@abundanceList);
		} else {
			@colorList = map {$_ > 1 ? color(1, 1 / $_, 0) : color($_, 1, 0)} @abundanceList;
		}
		foreach my $sampleIndex (0 .. $#sampleList) {
			my ($sample, $abundance, $color) = ($sampleList[$sampleIndex], $abundanceList[$sampleIndex], $colorList[$sampleIndex]);
			if(defined($_ = $sampleGeneTaxonAbundanceHash{$sample}->{$feature})) {
				my %taxonAbundanceHash = %$_;
				my @taxonAbundanceList = sort {$b->[1] <=> $a->[1] || $b->[0] cmp $a->[0]} map {[$_, $taxonAbundanceHash{$_}]} keys %taxonAbundanceHash;
				my $taxonAbundances = join('<br>', map {sprintf('%s %.3f', @$_)} @taxonAbundanceList);
				my @taxonFontSizeList = map {[substr($_->[0], 0, $taxonAbbreviationLength), $fontSize * ($_->[1] / $abundance)]} @taxonAbundanceList;
				my $taxons = join('', map {sprintf('<div style="font-size: %fpx; line-height: %fpx;">%s</div>', $_->[1], $_->[1], $_->[0])} @taxonFontSizeList);
				push(@tdList, "<td style=\"text-align: center; background-color: $color;\">$taxons<div class=\"popup\">$taxonAbundances</div></td>");
			} else {
				push(@tdList, "<td style=\"text-align: center; background-color: $color;\"></td>");
			}
		}
		print join('', "<tr>", @tdList, "</tr>"), "\n";
	}
	print "</table></td>\n";
	print "<td><table class=\"heatmap\">\n";
	if($scale) {
		foreach(map {$maximumScale * (5 - $_) / 5} 0 .. 4) {
			my $color = color(1, ($maximumScale - $_) / $maximumScale, 0);
			print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$_</td></tr>\n";
		}
		foreach(map {$maximumScale * (0 - $_) / 5} 0 .. 5) {
			my $color = color(($maximumScale - abs($_)) / $maximumScale, 1, 0);
			print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$_</td></tr>\n";
		}
	} else {
		foreach(map {10 - $_ * 2} 0 .. 4) {
			my $color = color(1, 1 / $_, 0);
			print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$_</td></tr>\n";
		}
		foreach(map {(5 - $_) / 5} 0 .. 5) {
			my $color = color($_, 1, 0);
			print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$_</td></tr>\n";
		}
	}
	print "</table></td>\n";
	print "</tr></table>\n";

	print <<EOF;
</form>
</body>
</html>
EOF
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

sub color {
	my ($red, $green, $blue) = @_;
	(my $color = sprintf('#%02x%02x%02x', int($red * 255), int($green * 255), int($blue * 255))) =~ tr/a-z/A-Z/;
	return $color;
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
