#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Bio::DB::Taxonomy;
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

my @refseqCategoryList = ('reference genome', 'representative genome', 'na');
my %refseqCategoryIndexHash = map {$refseqCategoryList[$_] => $_} 0 .. $#refseqCategoryList;

my @assemblyLevelList = ('Complete Genome', 'Chromosome', 'Scaffold', 'Contig');
my %assemblyLevelIndexHash = map {$assemblyLevelList[$_] => $_} 0 .. $#assemblyLevelList;

GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	't=s' => \(my $taxonIds = '2,2157,4751'),
	'u' => \(my $includeUnclassifiedTaxons = ''),
);
if($help) {
	die <<EOF;

Usage:   perl MetaPrism_data.taxon.pl [options]

Options: -h       display this help message
         -r       redownload data
         -t STR   NCBI taxonomy ID [$taxonIds]
         -u       include unclassified taxons

EOF
}
my @columnList = ();
{
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

		sub getUnclassifiedTaxonName {
			my ($taxonId) = @_;
			my $taxon = $db->get_taxon(-taxonid => $taxonId);
			do {
				return $taxon->scientific_name if($taxon->scientific_name =~ /^unclassified /);
			} while(defined($taxon = $taxon->ancestor));
			return '';
		}
	}
	open(my $writer, "| sort -t '\t' -k1,1n -k2,2n -k3,3n | cut -f4- > $dataPath/assembly_summary_refseq.sorted.txt");
	open(my $reader, "wget --no-verbose -O - ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt |");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^#/) {
			@columnList = split(/\t/, $line, -1) if($line =~ s/^# *//);
			next;
		}
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		if(defined($taxonIdHash{$tokenHash{'species_taxid'}}) && $tokenHash{'release_type'} eq 'Major' && $tokenHash{'genome_rep'} eq 'Full') {
			print $writer join("\t", $tokenHash{'species_taxid'}, $refseqCategoryIndexHash{$tokenHash{'refseq_category'}}, $assemblyLevelIndexHash{$tokenHash{'assembly_level'}}, @tokenHash{@columnList}), "\n";
		}
	}
	close($reader);
	close($writer);
}
{
	open(my $writer, "> $dataPath/species_assembly.fasta");
	open(my $reader, "$dataPath/assembly_summary_refseq.sorted.txt");
	my %topTokenHash = ();
	$topTokenHash{'species_taxid'} = '';
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		%topTokenHash = %tokenHash if($tokenHash{'species_taxid'} ne $topTokenHash{'species_taxid'});
		if(scalar(grep {$tokenHash{$_} ne $topTokenHash{$_}} 'refseq_category', 'assembly_level') == 0) {
			if($includeUnclassifiedTaxons eq '' && $tokenHash{'refseq_category'} eq 'na' && (my $unclassifiedTaxonName = getUnclassifiedTaxonName($tokenHash{'species_taxid'})) ne '') {
				print STDERR "Exclude $unclassifiedTaxonName\n";
				next;
			}
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
	}
	close($reader);
	close($writer);
#	system("rm $dataPath/assembly_summary_refseq.sorted.txt");
}

system("minimap2 -d $dataPath/species_assembly.mni -x asm5 $dataPath/species_assembly.fasta");
