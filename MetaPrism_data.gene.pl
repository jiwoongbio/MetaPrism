#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Bio::DB::Taxonomy;
use IPC::Open2;
use Time::HiRes;
use URI::Escape;
use XML::LibXML;
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

my @dataURLList = (
	'https://cdc.biohpc.swmed.edu/FMAP/FMAP_data',
	'https://qbrc.swmed.edu/FMAP/FMAP_data',
);

# NCBI E-utility
my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $apiKey = $ENV{'NCBI_API_KEY'};

GetOptions(
	'h' => \(my $help = ''),
	'p=i' => \(my $threads = 1),
	'b' => \(my $databaseOrthology = ''),
	'a' => \(my $databaseAMR = ''),
	'r' => \(my $redownload = ''),
	't=s' => \(my $taxonIds = '2,2157,4751'),
	'i=s' => \(my $unirefIdentity = '90'),
);
if($help) {
	die <<EOF;

Usage:   perl MetaPrism_data.gene.pl [options]

Options: -h       display this help message
         -p INT   number of threads [$threads]
         -r       redownload data
         -b       build orthology database
         -a       build AMR gene database
         -t STR   NCBI taxonomy ID [$taxonIds]
         -i INT   UniRef identity [$unirefIdentity]

EOF
}
my $database = '';
if($databaseOrthology) {
	my @taxonIdList = ();
	if($taxonIds ne '') {
		@taxonIdList = sort {$a <=> $b} eval($taxonIds);
		@taxonIdList = @taxonIdList[0, grep {$taxonIdList[$_ - 1] != $taxonIdList[$_]} 1 .. $#taxonIdList];
	}
	my %taxonIdHash = ();
	if(@taxonIdList) {
		if(not -r "$dataPath/nodes.dmp" or not -r "$dataPath/names.dmp" or $redownload) {
			my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
			my $file = "$dataPath/taxdump.tar.gz";
			system("wget --no-verbose --no-check-certificate -O $file $URL") if(not -r $file or $redownload);
			system("cd $dataPath; tar -zxf taxdump.tar.gz nodes.dmp");
			system("cd $dataPath; tar -zxf taxdump.tar.gz names.dmp");
			system("rm -f $dataPath/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
		}
		my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataPath, -nodesfile => "$dataPath/nodes.dmp", -namesfile => "$dataPath/names.dmp");
		my @taxonList = map {$db->get_taxon(-taxonid => $_)} @taxonIdList;
		setTaxonIdHash($_) foreach(@taxonList);

		sub setTaxonIdHash {
			my ($taxon) = @_;
			return if(defined($taxonIdHash{$taxon->id}));
			$taxonIdHash{$taxon->id} = $taxon;
			my @taxonList = $db->each_Descendent($taxon);
			setTaxonIdHash($_) foreach(@taxonList);
		}
	}
	$database = join('_', 'orthology', "uniref$unirefIdentity", @taxonIdList);
	$database = join('.', $database, getTimeString());
	my %unirefOrthologyCountHash = ();
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2 | uniq | cut -f2-");
	open(my $writerLog, "> $dataPath/$database.log");
	{
		my $URL = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz';
		my $file = "$dataPath/idmapping.dat.gz";
		system("wget --no-verbose --no-check-certificate -O $file $URL") if(not -r $file or $redownload);
		my %tokenHash = ('UniProtKB-AC' => '');
		open(my $reader, "gzip -dc $file | sort -t '\t' -k1,1 |");
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
			if($tokenList[0] ne $tokenHash{'UniProtKB-AC'}) {
				printGeneUniref() if($tokenHash{'UniProtKB-AC'} ne '');
				%tokenHash = ('UniProtKB-AC' => $tokenList[0]);
			}
			push(@{$tokenHash{$tokenList[1]}}, $tokenList[2]);
		}
		printGeneUniref() if($tokenHash{'UniProtKB-AC'} ne '');
		close($reader);

		sub printGeneUniref {
			my ($taxonIdList, $geneList, $unirefList) = @tokenHash{'NCBI_TaxID', 'KEGG', "UniRef$unirefIdentity"};
			if(scalar(keys %taxonIdHash) == 0 || (defined($taxonIdList) && grep {defined($taxonIdHash{$_})} @$taxonIdList)) {
				if(defined($unirefList) && defined($geneList)) {
					print $writer join("\t", $_->[1], $_->[0], join(',', @$unirefList)), "\n" foreach(map {[$_, split(/:/, $_)]} @$geneList);
				}
			}
		}
	}
	close($writer);
	{
		my $orgCode = '';
		my %geneUnirefHash = ();
		while(my $line = <$reader>) {
			chomp($line);
			my ($gene, $uniref) = split(/\t/, $line, -1);
			if($gene !~ /^$orgCode:/) {
				addUnirefOrthology() if($orgCode ne '');
				($orgCode = $gene) =~ s/:.*$//;
				%geneUnirefHash = ();
			}
			$geneUnirefHash{$gene}->{$uniref} = 1;
		}
		addUnirefOrthology() if($orgCode ne '');
		close($reader);

		sub addUnirefOrthology {
			my %geneOrthologyHash = ();
			open(my $reader, "wget --no-verbose --no-check-certificate -O - http://rest.kegg.jp/link/ko/$orgCode |");
			while(my $line = <$reader>) {
				chomp($line);
				my ($gene, $orthology) = split(/\t/, $line, -1);
				$orthology =~ s/^ko://;
				$geneOrthologyHash{$gene}->{$orthology} = 1;
			}
			close($reader);
			foreach my $gene (keys %geneUnirefHash) {
				my @unirefList = sort keys %{$geneUnirefHash{$gene}};
				my @orthologyList = ();
				@orthologyList = sort keys %$_ if(defined($_ = $geneOrthologyHash{$gene}));
				if(scalar(@unirefList) == 1 && scalar(@orthologyList) == 1) {
					$unirefOrthologyCountHash{$unirefList[0]}->{$orthologyList[0]} += 1;
				} else {
					print $writerLog join("\t", $gene, join(',', @unirefList), join(',', @orthologyList)), "\n";
				}
			}
		}
	}
	{
		my $URL = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref$unirefIdentity/uniref$unirefIdentity.fasta.gz";
		my $file = "$dataPath/uniref$unirefIdentity.fasta.gz";
		system("wget --no-verbose --no-check-certificate -O $file $URL") if(not -r $file or $redownload);
		open(my $reader, "gzip -dc $file |");
		open(my $writer, "> $dataPath/$database.fasta");
		my $printSequence = '';
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ /^>(\S*)/) {
				my $uniref = $1;
				$printSequence = '';
				if(defined(my $orthologyCountHash = $unirefOrthologyCountHash{$uniref})) {
					if(scalar(my @orthologyList = sort keys %$orthologyCountHash) == 1) {
						print $writer ">$orthologyList[0]_$uniref\n";
						$printSequence = 1;
					} else {
						print $writerLog join("\t", $uniref, $_, $orthologyCountHash->{$_}), "\n" foreach(@orthologyList);
					}
				}
			} elsif($printSequence) {
				print $writer "$line\n";
			}
		}
		close($reader);
		close($writer);
	}
	close($reader);
	waitpid($pid, 0);
	close($writerLog);
	{
		open(my $reader, 'wget --no-verbose --no-check-certificate -O - http://rest.kegg.jp/list/ko |');
		open(my $writer, "> $dataPath/$database.definition.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my ($orthology, $definition) = split(/\t/, $line, -1);
			$orthology =~ s/^ko://;
			print $writer join("\t", $orthology, $definition), "\n";
		}
		close($reader);
		close($writer);
	}
	{
		open(my $writer, "> $dataPath/database");
		print $writer "$database\n";
		close($writer);
	}
} elsif($databaseAMR) {
	$database = 'AMR';
	my %geneDefinitionCountHash = ();
	{
		open(my $reader, "wget --no-verbose --no-check-certificate -O - ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt |");
		open(my $writer, "| sort -t '\t' -k1,1 -k2,2 | uniq > $dataPath/$database.txt");
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			if($tokenHash{'type'} eq 'AMR' && ((my $protein = $tokenHash{'refseq_protein_accession'}) ne '') && ((my $gene = $tokenHash{'gene_family'}) ne '')) {
				print $writer join("\t", $protein, $gene), "\n";
				if((my $definition = $tokenHash{'product_name'}) ne '') {
					$geneDefinitionCountHash{$gene}->{$definition} += 1;
				}
			}
		}
		close($reader);
		close($writer);
	}
	{
		open(my $reader, "$dataPath/$database.txt");
		open(my $writer, "> $dataPath/$database.fasta");
		while(my $line = <$reader>) {
			chomp($line);
			my ($protein, $gene) = split(/\t/, $line, -1);
			my @idList = esearch('protein', sprintf('%s[Accession]', $protein));
			print $writer efetch('protein', $_, 'fasta', 'text') foreach(@idList);
		}
		close($reader);
		close($writer);
	}
	{
		open(my $writer, "> $dataPath/$database.definition.txt");
		foreach my $gene (sort keys %geneDefinitionCountHash) {
			my %definitionCountHash = %{$geneDefinitionCountHash{$gene}};
			my @definitionList = sort {$definitionCountHash{$b} <=> $definitionCountHash{$a} || $a cmp $b} keys %definitionCountHash;
			print $writer join("\t", $gene, $definitionList[0]), "\n";
		}
		close($writer);
	}
	{
		open(my $writer, "> $dataPath/database");
		print $writer "$database\n";
		close($writer);
	}
} else {
	downloadDataFile('database');
	open(my $reader, "$dataPath/database") or die "Can't open '$dataPath/database': $!";
	chomp($database = <$reader>);
	close($reader);
	downloadDataFile("$database.fasta.gz");
	system("gzip -df $dataPath/$database.fasta.gz");
	downloadDataFile("$database.definition.txt");
}
system("diamond makedb --threads $threads --in $dataPath/$database.fasta --db $dataPath/$database.dmnd");

sub getTimeString {
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
	return sprintf('%04d%02d%02d%02d%02d%02d', $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
}

sub downloadDataFile {
	foreach my $file (@_) {
		foreach my $dataURL (@dataURLList) {
			system("wget --no-verbose --no-check-certificate -O $dataPath/$file $dataURL/$file") if(not -r "$dataPath/$file" or $redownload);
		}
	}
}

sub efetch {
	my ($db, $id, $rettype, $retmode) = @_;
	my $output = '';
	while($output eq '') {
		if(defined($apiKey)) {
			$output = `wget --no-verbose --no-check-certificate -O - '$baseURL/efetch.fcgi?db=$db&id=$id&rettype=$rettype&retmode=$retmode&api_key=$apiKey'`;
			Time::HiRes::usleep(105);
		} else {
			$output = `wget --no-verbose --no-check-certificate -O - '$baseURL/efetch.fcgi?db=$db&id=$id&rettype=$rettype&retmode=$retmode'`;
			Time::HiRes::usleep(350);
		}
	}
	return $output;
}

sub esearch {
	my ($db, $term, $mindate, $maxdate) = @_;
	my $encodedTerm = uri_escape($term);
	$mindate = '' unless(defined($mindate));
	$maxdate = '' unless(defined($maxdate));
	my @idList = ();
	my $retmax = 500;
	for(my $retstart = 0, my ($count, $web, $key) = ($retmax); $retstart < $count;) {
		my $xmlString = '';
		unless(defined($web) && defined($key)) {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose --no-check-certificate -O - '$baseURL/esearch.fcgi?db=$db&term=$encodedTerm&usehistory=y&retmax=$retmax&retstart=$retstart&mindate=$mindate&maxdate=$maxdate&api_key=$apiKey'`;
			} else {
				$xmlString = `wget --no-verbose --no-check-certificate -O - '$baseURL/esearch.fcgi?db=$db&term=$encodedTerm&usehistory=y&retmax=$retmax&retstart=$retstart&mindate=$mindate&maxdate=$maxdate'`;
			}
		} else {
			if(defined($apiKey)) {
				$xmlString = `wget --no-verbose --no-check-certificate -O - '$baseURL/esearch.fcgi?WebEnv=$web&query_key=$key&retmax=$retmax&retstart=$retstart&api_key=$apiKey'`;
			} else {
				$xmlString = `wget --no-verbose --no-check-certificate -O - '$baseURL/esearch.fcgi?WebEnv=$web&query_key=$key&retmax=$retmax&retstart=$retstart'`;
			}
		}
		if($xmlString =~ /<\/eSearchResult>\n?$/) {
			my $dom = XML::LibXML->load_xml(string => $xmlString);
			my $root = $dom->documentElement();
			($count) = map {$_->textContent} getChildNodeList($root, 'Count');
			unless(defined($web) && defined($key)) {
				($web) = map {$_->textContent} getChildNodeList($root, 'WebEnv');
				($key) = map {$_->textContent} getChildNodeList($root, 'QueryKey');
			}
			push(@idList, map {$_->textContent} getChildNodeList($root, 'IdList', 'Id'));
			$retstart += $retmax;
		}
		if(defined($apiKey)) {
			Time::HiRes::usleep(105);
		} else {
			Time::HiRes::usleep(350);
		}
	}
	return @idList;
}

sub getChildNodeList {
	my ($node, @childNodeTagNameList) = @_;
	my @childNodeList = ();
	if(@childNodeTagNameList) {
		foreach my $childNode ($node->getChildrenByTagName(shift @childNodeTagNameList)) {
			push(@childNodeList, getChildNodeList($childNode, @childNodeTagNameList));
		}
	} else {
		push(@childNodeList, $node);
	}
	return @childNodeList;
}
