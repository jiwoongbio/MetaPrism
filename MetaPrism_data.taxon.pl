#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Bio::DB::Taxonomy;
use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

my @assemblySummaryFilePathList = (
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt',
);

my @refseqCategoryList = ('reference genome', 'representative genome', 'na');
my %refseqCategoryIndexHash = map {$refseqCategoryList[$_] => $_} 0 .. $#refseqCategoryList;

my @assemblyLevelList = ('Complete Genome', 'Chromosome', 'Scaffold', 'Contig');
my %assemblyLevelIndexHash = map {$assemblyLevelList[$_] => $_} 0 .. $#assemblyLevelList;

GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	't=s' => \(my $taxonIds = '2,2157,4751'),
);
if($help) {
	die <<EOF;

Usage:   perl MetaPrism_data.taxon.pl [options]

Options: -h       display this help message
         -r       redownload data
         -t STR   NCBI taxonomy ID [$taxonIds]

EOF
}
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
		system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
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

	sub isUnclassifiedOrUncultured {
		my ($taxonId) = @_;
		my $taxon = $db->get_taxon(-taxonid => $taxonId);
		while(defined($taxon)) {
			return 1 if(defined($taxon->scientific_name) && ($taxon->scientific_name =~ /^unclassified / || $taxon->scientific_name =~ /^uncultured / || $taxon->scientific_name =~ /^Candidatus / || $taxon->scientific_name =~ / \S*CAG\S*$/));
			$taxon = $taxon->ancestor;
		}
		return '';
	}
}
{
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1n -k2,2n -k3,3n -k4,4n -k5,5nr | cut -f6-");
	my @headerLineList = ();
	my @outputColumnList = ();
	foreach my $assemblySummaryFilePathIndex (0 .. $#assemblySummaryFilePathList) {
		my $assemblySummaryFilePath = $assemblySummaryFilePathList[$assemblySummaryFilePathIndex];
		open(my $reader, "wget --no-verbose -O - $assemblySummaryFilePath |");
		my $line;
		chomp($line = <$reader>);
		push(@headerLineList, $line) if($assemblySummaryFilePathIndex == 0);
		chomp($line = <$reader>);
		push(@headerLineList, $line) if($assemblySummaryFilePathIndex == 0);
		$line =~ s/^# ?//;
		my @columnList = split(/\t/, $line, -1);
		@outputColumnList = @columnList if($assemblySummaryFilePathIndex == 0);
		while($line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			next unless(defined($taxonIdHash{$tokenHash{'species_taxid'}}));
			next if($tokenHash{'ftp_path'} eq 'na');
			my $refseqCategoryIndex = $refseqCategoryIndexHash{$tokenHash{'refseq_category'}};
			my $assemblyLevelIndex = $assemblyLevelIndexHash{$tokenHash{'assembly_level'}};
			(my $assemblyAccessionNumber = $tokenHash{'assembly_accession'}) =~ s/^[^0-9]*//;
			my $assemblyAccessionNumberVersion = 0;
			$assemblyAccessionNumberVersion = $1 if($assemblyAccessionNumber =~ s/\.([0-9]+)$//);
			print $writer join("\t", $refseqCategoryIndex, $assemblyLevelIndex, $assemblySummaryFilePathIndex, $assemblyAccessionNumber, $assemblyAccessionNumberVersion, @tokenHash{@outputColumnList}), "\n";
		}
		close($reader);
	}
	close($writer);
	{
		open(my $writer, "> $dataPath/assembly_summary.sorted.txt");
		print $writer $_, "\n" foreach(@headerLineList);
		my %taxonIdHash = ();
		my %assemblyAccessionHash = ();
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@outputColumnList} = split(/\t/, $line, -1);
			next if($tokenHash{'refseq_category'} eq 'na' && defined($taxonIdHash{$tokenHash{'species_taxid'}}));
			$taxonIdHash{$tokenHash{'species_taxid'}} = 1;
			next if($assemblyAccessionHash{$tokenHash{'assembly_accession'}});
			next if(isUnclassifiedOrUncultured($tokenHash{'species_taxid'}));
			print $writer join("\t", @tokenHash{@outputColumnList}), "\n";
			$assemblyAccessionHash{$tokenHash{'assembly_accession'}} = 1;
			$assemblyAccessionHash{$tokenHash{'gbrs_paired_asm'}} = 1 if($tokenHash{'paired_asm_comp'} eq 'identical');
		}
		close($writer);
	}
	close($reader);
	waitpid($pid, 0);
}
{
	open(my $writer, "> $dataPath/species_assembly.fasta");
	open(my $reader, "$dataPath/assembly_summary.sorted.txt");
	my $line;
	chomp($line = <$reader>);
	chomp($line = <$reader>);
	$line =~ s/^# ?//;
	my @columnList = split(/\t/, $line, -1);
	while($line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		my $ftpPath = $tokenHash{'ftp_path'};
		(my $file = "$ftpPath\_genomic.fna.gz") =~ s/^.*\///;
		open(my $reader, "wget --no-verbose -O - $ftpPath/$file | gzip -d |");
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ s/^>//) {
				$line = ">$tokenHash{'species_taxid'}|$tokenHash{'assembly_accession'}|$line";
			}
			print $writer "$line\n";
		}
		close($reader);
	}
	close($reader);
	close($writer);
}

system("minimap2 -d $dataPath/species_assembly.mni -x asm5 $dataPath/species_assembly.fasta");
