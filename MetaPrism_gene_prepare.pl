# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/MetaPrism_data";
system("mkdir -p $dataPath");

my $dataURL = 'https://qbrc.swmed.edu/FMAP/FMAP_data';

GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'k' => \(my $downloadPrebuiltKEGG = ''),
	'a' => \(my $ardbDatabase = ''),
	'b' => \(my $betalactamasesDatabase = ''),
);
if($help) {
	die <<EOF;

Usage:   perl MetaPrism_gene_preapare.pl [options]

Options: -h       display this help message
         -r       redownload data
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -k       download prebuilt KEGG files
         -a       download ARDB database
         -b       download beta-lactamase database

EOF
}
die "ERROR in $0: 'wget' is not executable.\n" unless(-x getCommandPath('wget'));
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR in $0: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper =~ /^diamond/ || $mapper =~ /^usearch/);
die "ERROR in $0: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

# Database
my ($database, $databasePath) = ('', '');
if(-r "$dataPath/database") {
	$database = loadDefaultDatabase();
	$databasePath = "$dataPath/$database";
} else {
	if($ardbDatabase) {
		downloadFile('ARDB.database');
		system("mv $dataPath/ARDB.database $dataPath/database");
	} elsif($betalactamasesDatabase) {
		downloadFile('betalactamases.database');
		system("mv $dataPath/betalactamases.database $dataPath/database");
	} else {
		downloadFile('database');
	}
	$database = loadDefaultDatabase();
	downloadFile("$database.fasta.gz");
	system("gzip -df $dataPath/$database.fasta.gz");
	$databasePath = "$dataPath/$database";
	if($ardbDatabase || $betalactamasesDatabase) {
		downloadFile("$database.txt");
		downloadFile("$database.definition.txt");
	}
}
if(-r "$databasePath.fasta") {
	generateSequenceLengthFile("$databasePath.length.txt", "$databasePath.fasta");
	if($mapper =~ /^diamond/) {
		system("$mapperPath makedb --in $databasePath.fasta --db $databasePath.dmnd 1>&2");
	}
	if($mapper =~ /^usearch/) {
		system("$mapperPath -makeudb_ublast $databasePath.fasta -output $databasePath.udb 1>&2");
	}
} else {
	die "ERROR in $0: The database is not available.\n";
}

# KEGG
if($database =~ /^orthology_uniref/) {
	my $file = 'KEGG_orthology.txt';
	if($downloadPrebuiltKEGG) {
		downloadFile($file);
	} elsif(not -r "$dataPath/$file" or $redownload) {
		open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/ko |');
		open(my $writer, "> $dataPath/$file");
		while(my $line = <$reader>) {
			chomp($line);
			my ($orthology, $definition) = split(/\t/, $line);
			$orthology =~ s/^ko://;
			print $writer join("\t", $orthology, $definition), "\n";
		}
		close($reader);
		close($writer);
	}
}

sub downloadFile {
	foreach my $file (@_) {
		system("wget --no-verbose -O $dataPath/$file $dataURL/$file") if(not -r "$dataPath/$file" or $redownload);
	}
}

sub loadDefaultDatabase {
	open(my $reader, "$dataPath/database");
	chomp(my $line = <$reader>);
	close($reader);
	return $1 if($line =~ /^([A-Za-z0-9_.]+)/);
}

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}

sub generateSequenceLengthFile {
	my ($sequenceLengthFile, $fastaFile) = @_;
	open(my $writer, "> $sequenceLengthFile");
	my ($sequenceName, $sequenceLength) = ('', 0);
	open(my $reader, $fastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(.*)$/) {
			(my $nextSequenceName = $1) =~ s/ .*$//;
			print $writer join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
			($sequenceName, $sequenceLength) = ($nextSequenceName, 0);
		} else {
			$sequenceLength += length($line);
		}
	}
	close($reader);
	print $writer join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
	close($writer);
}
