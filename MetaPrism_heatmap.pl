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

my @geneFileList = ();
GetOptions(
	'h' => \(my $help = ''),
	'A=s' => \(my $abundanceColumn = 'meanDepth/genome'),
	'R=s' => \(my $taxonRank = 'genus'),
	'g=s' => \@geneFileList,
	't=i' => \(my $taxonAbbreviationLength = 4),
	'f=i' => \(my $fontSize = 15),
	'w=i' => \(my $tdWidth = 60),
	'S' => \(my $disableClusteringSample = ''),
	'G' => \(my $disableClusteringGene = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_heatmap.pl [options] [sample=]abundance.txt [...] > heatmap.html

Options: -h       display this help message
         -A STR   abundance column [$abundanceColumn]
         -R STR   taxon rank [$taxonRank]
         -g FILE  gene file
         -t INT   taxon abbreviation length [$taxonAbbreviationLength]
         -f INT   HTML font size [$fontSize]
         -w INT   HTML table cell width [$tdWidth]
         -S       disable clustering sample
         -G       disable clustering gene

EOF
}

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

my %sampleGeneAbundanceHash = ();
my %sampleGeneTaxonAbundanceHash = ();
my %geneDefinitionHash = ();
@sampleFileList = map {[$_->[0], $_->[-1]]} map {[split(/=/, $_, 2)]} @sampleFileList;
foreach(@sampleFileList) {
	my ($sample, $file) = @$_;
	if(-s $file) {
		open(my $reader, $file);
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			$sampleGeneAbundanceHash{$sample}->{$tokenHash{'gene'}} += $tokenHash{$abundanceColumn};
			if((my $taxonName = getTaxonName($tokenHash{'taxon'})) ne '') {
				$sampleGeneTaxonAbundanceHash{$sample}->{$tokenHash{'gene'}}->{$taxonName} += $tokenHash{$abundanceColumn};
			}
			$geneDefinitionHash{$tokenHash{'gene'}} = $tokenHash{'definition'} if(defined($tokenHash{'definition'}) && $tokenHash{'definition'} ne '');
		}
		close($reader);
	}
}
my @sampleList = sort keys %sampleGeneAbundanceHash;
my @geneList = ();
{
	my %geneHash = ();
	foreach my $sample (@sampleList) {
		$geneHash{$_} = 1 foreach(keys %{$sampleGeneAbundanceHash{$sample}});
	}
	@geneList = sort keys %geneHash;
}
if(@geneFileList) {
	my %geneHash = ();
	foreach my $geneFile (@geneFileList) {
		open(my $reader, $geneFile);
		while(my $line = <$reader>) {
			chomp($line);
			my ($gene) = split(/\t/, $line);
			$geneHash{$gene} = 1;
		}
		close($reader);
	}
	@geneList = grep {$geneHash{$_}} @geneList;
}

if(@geneList) {
	my $R = Statistics::R->new();
	$R->run('x <- data.frame()');
	foreach my $geneIndex (0 .. $#geneList) {
		foreach my $sampleIndex (0 .. $#sampleList) {
			if(defined(my $abundance = $sampleGeneAbundanceHash{$sampleList[$sampleIndex]}->{$geneList[$geneIndex]})) {
				$R->set(sprintf('x[%d, %d]', $sampleIndex + 1, $geneIndex + 1), $abundance);
			} else {
				$R->set(sprintf('x[%d, %d]', $sampleIndex + 1, $geneIndex + 1), 0);
			}
		}
	}
	foreach my $sampleIndex (0 .. $#sampleList) {
		$R->set(sprintf('rownames(x)[%d]', $sampleIndex + 1), $sampleList[$sampleIndex]);
	}
	foreach my $geneIndex (0 .. $#geneList) {
		$R->set(sprintf('colnames(x)[%d]', $geneIndex + 1), $geneList[$geneIndex]);
	}
	if($disableClusteringSample eq '') {
		$R->run('hc.sample <- hclust(as.dist(1 - cor(t(x), method = "spearman")))');
		@sampleList = @{$R->get('hc.sample$labels[hc.sample$order]')};
	}
	if($disableClusteringGene eq '') {
		$R->run('hc.gene <- hclust(as.dist(1 - cor(x, method = "spearman")))');
		@geneList = @{$R->get('hc.gene$labels[hc.gene$order]')};
	}
	$R->stop();
}

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
	my @tdList = ("<td></td>");
	foreach my $sample (@sampleList) {
		push(@tdList, "<td style=\"text-align: center; width: ${tdWidth}px;\">$sample</td>");
	}
	print join('', "<tr>", @tdList, "</tr>"), "\n";
}
foreach my $gene (@geneList) {
	my @tdList = ();
	if(defined(my $definition = $geneDefinitionHash{$gene})) {
		push(@tdList, "<td>$gene<div class=\"popup\">$definition</div></td>");
	} else {
		push(@tdList, "<td>$gene</td>");
	}
	foreach my $sample (@sampleList) {
		my $abundance = defined($_ = $sampleGeneAbundanceHash{$sample}->{$gene}) ? $_ : 0;
		my $color = '#FFFFFF';
		if($abundance > 1) {
			$color = color(1, 1 / $abundance, 0);
		} else {
			$color = color($abundance, 1, 0);
		}
		if(defined($_ = $sampleGeneTaxonAbundanceHash{$sample}->{$gene})) {
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
for(my $depth = 10; $depth > 0; $depth -= 2) {
	my $color = color(1, 1 / $depth, 0);
	print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$depth</td></tr>\n";
}
for(my $depth = 1; $depth >= 0; $depth = sprintf('%.1f', $depth - 0.2)) {
	my $color = color($depth, 1, 0);
	print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$depth</td></tr>\n";
}
print "</table></td>\n";
print "</tr></table>\n";

print <<EOF;
</form>
</body>
</html>
EOF

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
