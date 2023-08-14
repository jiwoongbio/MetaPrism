#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use List::Util qw(sum max min);
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 100),
	'P' => \(my $inputIsProteinFastaFile = ''),
	'C=s' => \@codonList,
	'S=s' => \(my $startCodons = 'GTG,ATG,CTG,TTG,ATA,ATC,ATT'),
	'T=s' => \(my $terminationCodons = 'TAG,TAA,TGA'),
	'l=i' => \(my $minimumTranslationLength = 10),
	'e=f' => \(my $evalue = 10),
	'k=i' => \(my $maximumTarget = 100),
	'c=f' => \(my $minimumCoverage = 0.8),
	'd=s' => \(my $diamond = 'diamond'),
	'D=s' => \(my $diamondDatabaseFile),
	'G=s' => \(my $geneDefinitionFile),
	'T=s' => \(my $proteinGeneFile),
	'm=s' => \(my $minimap2File),
	'r=s' => \(my $rank = 'species'),
	'taxonFile=s' => \(my $taxonFile = ''),
	'pafFile=s' => \(my $pafFile = ''),
	'printAllProteinMappings' => \(my $printAllProteinMappings = ''),
	'numberPerProteinMapping=i' => \(my $numberPerProteinMapping = 10000000),
	'fastaLineLength=i' => \(my $fastaLineLength = 80),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism.pl [options] genome.fasta > MetaPrism.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -P       input is protein fasta file
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [$startCodons]
         -T STR   comma-separated termination codons [$terminationCodons]
         -l INT   minimum translation length [$minimumTranslationLength]
         -e FLOAT maximum e-value to report alignments [$evalue]
         -k INT   maximum target [$maximumTarget]
         -c FLOAT minimum coverage [$minimumCoverage]
         -d FILE  diamond path [$diamond]
         -D FILE  diamond database file
         -G FILE  gene definition file
         -T FILE  protein-to-gene file
         -m FILE  minimap2 file
         -r STR   taxonomic rank [$rank]

EOF
}
my $temporaryPrefix = "$temporaryDirectory/$hostname.$$";
system("rm -f $temporaryPrefix.*");
{
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryPrefix.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryPrefix.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryPrefix.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
{
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. int(length($sequence) / 3) - 1);
	}
}
my %startCodonHash = map {$_ => 1} split(/,/, $startCodons);
my %terminationCodonHash = map {$_ => 1} split(/,/, $terminationCodons);
my $minimumLength = $minimumTranslationLength * 3;
unless(defined($diamondDatabaseFile)) {
	open(my $reader, "$dataPath/database") or die "Can't open '$dataPath/database': $!";
	chomp(my $database = <$reader>);
	close($reader);
	$diamondDatabaseFile = "$dataPath/$database.dmnd";
	unless(defined($geneDefinitionFile)) {
		$geneDefinitionFile = "$dataPath/$database.definition.txt" if(-r "$dataPath/$database.definition.txt");
	}
	unless(defined($proteinGeneFile)) {
		$proteinGeneFile = "$dataPath/$database.txt" if(-r "$dataPath/$database.txt");
	}
	unless(defined($minimap2File)) {
		$minimap2File = "$dataPath/species_assembly.mni" if(-r "$dataPath/species_assembly.mni");
	}
}
my %geneDefinitionHash = ();
if(defined($geneDefinitionFile)) {
	open(my $reader, ($geneDefinitionFile =~ /\.gz$/ ? "gzip -dc $geneDefinitionFile |" : $geneDefinitionFile)) or die "Can't open '$geneDefinitionFile': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my ($gene, $definition) = split(/\t/, $line, -1);
		$geneDefinitionHash{$gene} = $definition;
	}
	close($reader);
}
my %proteinGeneHash = ();
my %proteinMinimumPercentIdentityHash = ();
if(defined($proteinGeneFile)) {
	open(my $reader, ($proteinGeneFile =~ /\.gz$/ ? "gzip -dc $proteinGeneFile |" : $proteinGeneFile)) or die "Can't open '$proteinGeneFile': $!";
	while(my $line = <$reader>) {
		chomp($line);
		my ($protein, $gene, $minimumPercentIdentity) = split(/\t/, $line, -1);
		$proteinGeneHash{$protein} = $gene;
		$proteinMinimumPercentIdentityHash{$protein} = $minimumPercentIdentity;
	}
	close($reader);
}
my @columnList = $inputIsProteinFastaFile ? ('input') : ('chromosome', 'start', 'end', 'strand');
if(%geneDefinitionHash) {
	push(@columnList, 'gene', 'definition', 'protein', 'variant');
} else {
	push(@columnList, 'gene', 'protein', 'variant');
}
my @samMandatoryColumnList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
my @pafMandatoryColumnList = ('query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name', 'target_length', 'target_start', 'target_end', 'match', 'alignment_length', 'mapping_quality');
my ($genomeFastaFile) = @ARGV;
if($inputIsProteinFastaFile eq '') {
	open(my $writer, "> $temporaryPrefix.fasta");
	forkPrintParentWriter($writer);
	{
		open(my $reader, ($genomeFastaFile =~ /\.gz$/ ? "gzip -dc $genomeFastaFile |" : $genomeFastaFile));
		my ($chromosome, $sequence) = ('', '');
		my @chromosomeSequenceList = ();
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ /^>(\S*)/) {
				push(@chromosomeSequenceList, [$chromosome, $sequence]) if($sequence ne '');
				if(scalar(@chromosomeSequenceList) >= $numberPerThread) {
					if($threads == 1) {
						printChromosomeTranslationSequences(@chromosomeSequenceList);
					} else {
						forkPrintSubroutine(\&printChromosomeTranslationSequences, @chromosomeSequenceList);
					}
					@chromosomeSequenceList = ();
				}
				($chromosome, $sequence) = ($1, '');
				$chromosome =~ s/\|/_/g;
			} else {
				$sequence .= $line;
			}
		}
		push(@chromosomeSequenceList, [$chromosome, $sequence]) if($sequence ne '');
		if(@chromosomeSequenceList) {
			if($threads == 1) {
				printChromosomeTranslationSequences(@chromosomeSequenceList);
			} else {
				forkPrintSubroutine(\&printChromosomeTranslationSequences, @chromosomeSequenceList);
			}
		}
		close($reader);
		forkPrintWait();
	}
	forkPrintParentWriter();
	close($writer);

	sub printChromosomeTranslationSequences {
		foreach(@_) {
			my ($chromosome, $sequence) = @$_;
			my $sequenceLength = length($sequence);
			foreach my $frame (0 .. 2) {
				my @startIndexList = ();
				for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
					my $codon = substr($sequence, $index, 3);
					if($startCodonHash{$codon}) {
						push(@startIndexList, $index);
					} elsif(@startIndexList && $terminationCodonHash{$codon}) {
						printTranslationSequence($chromosome, $sequence, $sequenceLength, '+', $index + 3, @startIndexList);
						@startIndexList = ();
					}
				}
			}
			$sequence = getReverseComplementarySequence($sequence);
			foreach my $frame (0 .. 2) {
				my @startIndexList = ();
				for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
					my $codon = substr($sequence, $index, 3);
					if($startCodonHash{$codon}) {
						push(@startIndexList, $index);
					} elsif(@startIndexList && $terminationCodonHash{$codon}) {
						printTranslationSequence($chromosome, $sequence, $sequenceLength, '-', $index + 3, @startIndexList);
						@startIndexList = ();
					}
				}
			}
		}
	}

	sub printTranslationSequence {
		my ($chromosome, $sequence, $sequenceLength, $strand, $endIndex, @startIndexList) = @_;
		if((my $length = $endIndex - $startIndexList[0]) >= $minimumLength) {
			my ($start, $end) = ('', '');
			($start, $end) = ($startIndexList[0] + 1, $endIndex) if($strand eq '+');
			($start, $end) = (($sequenceLength - $endIndex) + 1, ($sequenceLength - $startIndexList[0])) if($strand eq '-');
			my @startList = map {($_ - $startIndexList[0]) / 3 + 1} @startIndexList;
			if((my $translationSequence = translate(substr($sequence, $startIndexList[0], $length))) =~ s/\*$//) {
				if(length($translationSequence) >= $minimumTranslationLength) {
					forkPrint('>', join('|', $chromosome, $start, $end, $strand, @startList), "\n");
					forkPrint("$translationSequence\n");
				}
			}
		}
	}
}
my %chromosomeHash = ();
{
	open(my $writer, ($inputIsProteinFastaFile ? "| sort -u > $temporaryPrefix.gene.txt" : "| sort -t '\t' -k1,1 -k2,2n -k3,3n -k4 | uniq > $temporaryPrefix.gene.txt"));
	open(my $reader, ($inputIsProteinFastaFile ? $genomeFastaFile : "$temporaryPrefix.fasta"));
	my $protein;
	my %proteinSequenceHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			if(scalar(keys %proteinSequenceHash) >= $numberPerProteinMapping) {
				printGenes();
				%proteinSequenceHash = ();
			}
			$protein = $1;
		} else {
			$proteinSequenceHash{$protein} .= $line;
		}
	}
	if(%proteinSequenceHash) {
		printGenes();
		%proteinSequenceHash = ();
	}
	close($reader);
	close($writer);
	system("rm $temporaryPrefix.fasta") if($inputIsProteinFastaFile eq '');

	sub printGenes {
		{
			open(my $writer, "| $diamond blastp --threads $threads --db $diamondDatabaseFile --out $temporaryPrefix.sam --outfmt 101 --evalue $evalue --unal 0 --tmpdir $temporaryDirectory --masking 0 --quiet --max-target-seqs $maximumTarget");
			foreach my $protein (keys %proteinSequenceHash) {
				print $writer ">$protein\n";
				writeFastaSequence($writer, $proteinSequenceHash{$protein});
			}
			close($writer);
		}
		{
			open(my $writer, "| sort -t '\t' -k4,4 -k1,1g -k2,2gr -k3,3gr | cut -f4- > $temporaryPrefix.sorted.sam");
			open(my $reader, "$temporaryPrefix.sam");
			while(my $line = <$reader>) {
				chomp($line);
				next if($line =~ /^@/);
				my %tokenHash = ();
				(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
				$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
				next if($tokenHash{'flag'} & 4);
				next if($tokenHash{'flag'} & 16);
				if(defined(my $minimumPercentIdentity = $proteinMinimumPercentIdentityHash{$tokenHash{'rname'}})) {
					next if($tokenHash{'ZI:i'} < $minimumPercentIdentity);
				}
				if(($tokenHash{'ZC:f'} = sum(0, $tokenHash{'MD:Z'} =~ /([0-9]+)/g) / $tokenHash{'ZL:i'}) >= $minimumCoverage) {
					push(@tagTypeValueList, "ZC:f:$tokenHash{'ZC:f'}");
					print $writer join("\t", @tokenHash{'ZE:f', 'AS:i', 'ZC:f'}, @tokenHash{@samMandatoryColumnList}, @tagTypeValueList), "\n";
				}
			}
			close($reader);
			close($writer);
			system("rm $temporaryPrefix.sam");
		}
		{
			open(my $reader, "$temporaryPrefix.sorted.sam");
			my @tokenHashList = ();
			while(my $line = <$reader>) {
				chomp($line);
				my %tokenHash = ();
				(@tokenHash{@samMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
				$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
				if(@tokenHashList && $tokenHash{'qname'} ne $tokenHashList[0]->{'qname'}) {
					printGene(@tokenHashList);
					@tokenHashList = ();
				}
				push(@tokenHashList, \%tokenHash);
			}
			if(@tokenHashList) {
				printGene(@tokenHashList);
				@tokenHashList = ();
			}
			close($reader);
			system("rm $temporaryPrefix.sorted.sam");
		}
	}

	sub printGene {
		my @tokenHashList = @_;
		my @queryEndList = ();
		foreach my $tokenHash (@tokenHashList) {
			my $querySequenceLength;
			if($inputIsProteinFastaFile) {
				$querySequenceLength = length($proteinSequenceHash{$tokenHash->{'qname'}});
			} else {
				my ($chromosome, $start, $end, $strand, @startList) = split(/\|/, $tokenHash->{'qname'});
				my $queryStartOffset;
				foreach my $queryStart (@startList) {
					last if($queryStart > $tokenHash->{'ZS:i'});
					$queryStartOffset = $queryStart - 1;
				}
				$start += $queryStartOffset * 3 if($strand eq '+');
				$end   -= $queryStartOffset * 3 if($strand eq '-');
				$tokenHash->{'qname'} = join('|', $chromosome, $start, $end, $strand);
				$tokenHash->{'ZS:i'} -= $queryStartOffset;
				$querySequenceLength = ($end - $start + 1) / 3 - 1;
			}
			my @variantList = getVariantList(@$tokenHash{'pos', 'cigar', 'MD:Z', 'seq', 'ZL:i', 'ZS:i'}, $querySequenceLength);
			$tokenHash->{'variant'} = \@variantList;
			my @queryStartEndList = ();
			$queryStartEndList[0]->[0] = 0;
			foreach my $index (0 .. $#variantList) {
				my ($queryPosition, $position, $referenceAA, $variantAA) = @{$variantList[$index]};
				$queryStartEndList[$index]->[1] = $queryPosition - 1;
				$queryStartEndList[$index + 1]->[0] = $queryPosition + ($variantAA =~ /^[A-Z]*$/ ? length($variantAA) : $variantAA);
			}
			$queryStartEndList[$#variantList + 1]->[1] = $querySequenceLength + 1;
			foreach(@queryStartEndList) {
				my ($queryStart, $queryEnd) = @$_;
				foreach my $queryPosition ($queryStart .. $queryEnd) {
					$queryEndList[$queryPosition] = $queryEnd if(!defined($queryEndList[$queryPosition]) || $queryEndList[$queryPosition] < $queryEnd);
				}
			}
		}
		foreach my $tokenHash (@tokenHashList) {
			next if($printAllProteinMappings eq '' && grep {$tokenHash->{$_} != $tokenHashList[0]->{$_}} 'ZE:f', 'AS:i', 'ZC:f');
			if($inputIsProteinFastaFile) {
				$tokenHash->{'input'} = $tokenHash->{'qname'};
			} else {
				@$tokenHash{'chromosome', 'start', 'end', 'strand'} = split(/\|/, $tokenHash->{'qname'});
				$chromosomeHash{$tokenHash->{'chromosome'}} = 1;
			}
			$tokenHash->{'protein'} = $tokenHash->{'rname'};
			if($tokenHash->{'protein'} =~ s/^(K[0-9]{5})_//) {
				$tokenHash->{'gene'} = $1;
			} elsif(defined(my $gene = $proteinGeneHash{$tokenHash->{'protein'}})) {
				$tokenHash->{'gene'} = $gene;
			} else {
				$tokenHash->{'gene'} = $tokenHash->{'protein'};
			}
			next if($tokenHash->{'gene'} eq 'off-target');
			$tokenHash->{'definition'} = $geneDefinitionHash{$tokenHash->{'gene'}};
			$tokenHash->{'definition'} = '' unless(defined($tokenHash->{'definition'}));
			my @variantList = ();
#			foreach my $variant (@{$tokenHash->{'variant'}}) {
#				my ($queryPosition, $position, $referenceAA, $variantAA) = @$variant;
#				if($variantAA eq '') {
#					push(@variantList, $variant) if(!defined($queryEndList[$queryPosition - 1]) || $queryEndList[$queryPosition - 1] < $queryPosition);
#				} else {
#					push(@variantList, $variant) if(!defined($queryEndList[$queryPosition]) || $queryEndList[$queryPosition] < $queryPosition + length($variantAA) - 1);
#				}
#			}
			@variantList = @{$tokenHash->{'variant'}};
			$tokenHash->{'variant'} = join(',', map {join('|', @$_[1, 2, 3])} @variantList);
			print $writer join("\t", @$tokenHash{@columnList}), "\n";
		}
	}

	sub getVariantList {
		my ($position, $cigar, $md, $sequence, $sequenceLength, $queryStart, $querySequenceLength) = @_;
		my @variantList = ();
		if($position - 1 > 0) {
			if($queryStart - 1 > 0) {
				push(@variantList, [1, 1, $position - 1, $queryStart - 1]);
			} else {
				push(@variantList, [1, 1, $position - 1, '']);
			}
		} else {
			if($queryStart - 1 > 0) {
				push(@variantList, [1, 1, '', $queryStart - 1]);
			}
		}
		my $index = 0;
		while($cigar =~ s/^([0-9]+)([MIDN])//) {
			my ($length, $operation) = ($1, $2);
			if($operation eq 'M') {
				if($md =~ s/^([0-9]+)// && $1 > 0) {
					if($length > $1) {
						$cigar = join('', $length - $1, 'M', $cigar);
						$length = $1;
					} elsif($1 > $length) {
						$md = join('', $1 - $length, $md);
					}
					$index += $length;
					$position += $length;
				} elsif($md =~ s/^([A-Z]+)0*//) {
					if($length > length($1)) {
						$cigar = join('', $length - length($1), 'M', $cigar);
						$length = length($1);
					} elsif(length($1) > $length) {
						$md = join('', substr($1, $length), $md);
					}
					push(@variantList, [$queryStart + $index, $position, $1, substr($sequence, $index, $length)]);
					$index += $length;
					$position += $length;
				} else {
					die "Failed to parse CIGAR/MD";
				}
			} elsif($operation eq 'I') {
					push(@variantList, [$queryStart + $index, $position, '', substr($sequence, $index, $length)]);
					$index += $length;
			} elsif($operation eq 'D') {
				if($md =~ s/^\^([A-Z]+)0*// && length($1) == $length) {
					push(@variantList, [$queryStart + $index, $position, $1, '']);
					$position += $length;
				} else {
					die "Failed to parse CIGAR/MD";
				}
			} elsif($operation eq 'N') {
				$position += $length;
			}
		}
		if($sequenceLength - ($position - 1) > 0) {
			if($querySequenceLength - ($queryStart + $index - 1) > 0) {
				push(@variantList, [$queryStart + $index, $position, $sequenceLength - ($position - 1), $querySequenceLength - ($queryStart + $index - 1)]);
			} else {
				push(@variantList, [$queryStart + $index, $position, $sequenceLength - ($position - 1), '']);
			}
		} else {
			if($querySequenceLength - ($queryStart + $index - 1) > 0) {
				push(@variantList, [$queryStart + $index, $position, '', $querySequenceLength - ($queryStart + $index - 1)]);
			}
		}
		die "Failed to parse CIGAR/MD" unless($cigar eq '' && $md eq '' && $index == length($sequence));
		return @variantList;
	}
}
if(%chromosomeHash && $minimap2File) {
	{
		open(my $writer, "> $temporaryPrefix.genome.fasta");
		open(my $reader, ($genomeFastaFile =~ /\.gz$/ ? "gzip -dc $genomeFastaFile |" : $genomeFastaFile));
		my ($chromosome, $sequence) = ('', '');
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ /^>(\S*)/) {
				if($sequence ne '' && $chromosomeHash{$chromosome}) {
					print $writer ">$chromosome\n";
					writeFastaSequence($writer, $sequence);
				}
				($chromosome, $sequence) = ($1, '');
				$chromosome =~ s/\|/_/g;
			} else {
				$sequence .= $line;
			}
		}
		if($sequence ne '' && $chromosomeHash{$chromosome}) {
			print $writer ">$chromosome\n";
			writeFastaSequence($writer, $sequence);
		}
		close($reader);
		close($writer);
	}
	{
		use Bio::DB::Fasta;
		my $db = Bio::DB::Fasta->new("$temporaryPrefix.genome.fasta");
		open(my $writer, "> $temporaryPrefix.fasta");
		{
			open(my $reader, "$temporaryPrefix.genome.fasta");
			print $writer $_ while(<$reader>);
			close($reader);
		}
		open(my $reader, "$temporaryPrefix.gene.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			print $writer ">$tokenHash{'chromosome'}|$tokenHash{'start'}|$tokenHash{'end'}\n";
			writeFastaSequence($writer, uc($db->seq($tokenHash{'chromosome'}, $tokenHash{'start'}, $tokenHash{'end'})));
		}
		close($reader);
		close($writer);
		system("rm $temporaryPrefix.genome.fasta");
		system("rm $temporaryPrefix.genome.fasta.index");
	}
	{
		open(my $writer, "| sort -t '\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5gr -k6,6gr | uniq > $temporaryPrefix.taxon.txt");
		open(my $reader, "minimap2 -t $threads --split-prefix $temporaryPrefix -x asm5 --secondary=no $minimap2File $temporaryPrefix.fasta |");
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^@/);
			my %tokenHash = ();
			(@tokenHash{@pafMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
#			$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
			(my $taxonId = $tokenHash{'target_name'}) =~ s/\|.*$//;
			my ($chromosome, $start, $end, $match, $identity) = ($tokenHash{'query_name'}, $tokenHash{'query_start'} + 1, $tokenHash{'query_end'}, $tokenHash{'match'}, $tokenHash{'match'} / $tokenHash{'alignment_length'});
			if($chromosome =~ s/\|([0-9]*)\|([0-9]*)$//) {
				($start, $end, $identity) = ($1, $2, $match / (($2 - $1) - ($end - $start) + $tokenHash{'alignment_length'}));
			}
			print $writer join("\t", $chromosome, $start, $end, $taxonId, $match, $identity), "\n";
		}
		close($reader);
		close($writer);
		system("rm $temporaryPrefix.fasta");
	}
	open(my $reader, "$temporaryPrefix.taxon.txt");
	{
		my @tokenList = getTokenList();
		my @tokenListList = ();
		sub getTokenList {
			while(my $line = <$reader>) {
				chomp($line);
				return split(/\t/, $line, -1);
			}
			return ();
		}
		sub getMatchIdentityTaxonIdList {
			my ($chromosome, $start, $end) = @_;
			@tokenListList = grep {$_->[0] eq $chromosome && $start <= $_->[2]} @tokenListList;
			while(@tokenList && $tokenList[0] lt $chromosome) {
				@tokenList = getTokenList;
			}
			while(@tokenList && $tokenList[0] eq $chromosome && $tokenList[1] <= $end) {
				push(@tokenListList, [@tokenList]) if($start <= $tokenList[2]);
				@tokenList = getTokenList;
			}
			my %matchIdentityTaxonIdHash = ();
			foreach(grep {$_->[1] <= $start && $end <= $_->[2]} @tokenListList) {
				my ($chromosome, $start, $end, $taxonId, $match, $identity) = @$_;
				$matchIdentityTaxonIdHash{$match}->{$identity}->{$taxonId} = 1;
			}
			if(%matchIdentityTaxonIdHash) {
				my $match = max(keys %matchIdentityTaxonIdHash);
				my $identity = max(keys %{$matchIdentityTaxonIdHash{$match}});
				my @taxonIdList = sort {$a <=> $b} keys %{$matchIdentityTaxonIdHash{$match}->{$identity}};
				@taxonIdList = @taxonIdList[0, grep {$taxonIdList[$_ - 1] != $taxonIdList[$_]} 1 .. $#taxonIdList];
				return ($match, $identity, @taxonIdList);
			} else {
				return ();
			}
		}
	}
	print join("\t", @columnList, 'taxonId', 'taxonName', 'taxonRank'), "\n";
	if(-s "$temporaryPrefix.gene.txt") {
		use Bio::DB::Taxonomy;
		my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataPath, -nodesfile => "$dataPath/nodes.dmp", -namesfile => "$dataPath/names.dmp");
		open(my $reader, "$temporaryPrefix.gene.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			my ($match, $identity, @taxonIdList) = getMatchIdentityTaxonIdList(@tokenHash{'chromosome', 'start', 'end'});
			my @rankTaxonIdList = getRankTaxonIdList($rank, @taxonIdList);
			if(scalar(@rankTaxonIdList) > 1) {
				my %taxonIdCountHash = ();
				my %taxonIdIndexHash = ();
				foreach my $taxonId (@rankTaxonIdList) {
					my $index = 0;
					my $taxon = $db->get_taxon(-taxonid => $taxonId);
					while(defined($taxon)) {
						$taxonIdCountHash{$taxon->id} += 1;
						$taxonIdIndexHash{$taxon->id} = $index if(!defined($taxonIdIndexHash{$taxon->id}) || $index < $taxonIdIndexHash{$taxon->id});
						$index += 1;
						$taxon = $taxon->ancestor;
					}
				}
				my %countTaxonIdListHash = ();
				push(@{$countTaxonIdListHash{$taxonIdCountHash{$_}}}, $_) foreach(keys %taxonIdCountHash);
				my %indexTaxonIdHash = ();
				$indexTaxonIdHash{$taxonIdIndexHash{$_}} = $_ foreach(@{$countTaxonIdListHash{max(keys %countTaxonIdListHash)}});
				my $taxonId = $indexTaxonIdHash{min(keys %indexTaxonIdHash)};
				my $taxon = $db->get_taxon(-taxonid => $taxonId);
				@tokenHash{'taxonId', 'taxonName', 'taxonRank'} = (join(',', @taxonIdList), $taxon->scientific_name, $taxon->rank);
			} elsif(scalar(@rankTaxonIdList) == 1) {
				my $taxonId = $rankTaxonIdList[0];
				my $taxon = $db->get_taxon(-taxonid => $taxonId);
				if(defined($taxon)) {
					@tokenHash{'taxonId', 'taxonName', 'taxonRank'} = (join(',', @taxonIdList), $taxon->scientific_name, $taxon->rank);
				} else {
					@tokenHash{'taxonId', 'taxonName', 'taxonRank'} = (join(',', @taxonIdList), '', '');
				}
			} else {
				@tokenHash{'taxonId', 'taxonName', 'taxonRank'} = ('', '', '');
			}
			print join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList, 'taxonId', 'taxonName', 'taxonRank'}), "\n";
		}
		close($reader);
		system("rm $temporaryPrefix.gene.txt");

		sub getRankTaxonIdList {
			my ($rank, @taxonIdList) = @_;
			my %taxonIdHash = ();
			foreach my $taxonId (@taxonIdList) {
				my $taxon = $db->get_taxon(-taxonid => $taxonId);
				while(defined($taxon)) {
					if($taxon->rank eq $rank) {
						$taxonIdHash{$taxon->id} = 1;
						last;
					}
					$taxon = $taxon->ancestor;
				}
			}
			return sort {$a <=> $b} keys %taxonIdHash;
		}
	}
	close($reader);
	if($taxonFile ne '') {
		system("mv $temporaryPrefix.taxon.txt $taxonFile");
	} else {
		system("rm $temporaryPrefix.taxon.txt");
	}
} else {
	print join("\t", @columnList), "\n";
	if(-s "$temporaryPrefix.gene.txt") {
		open(my $reader, "$temporaryPrefix.gene.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			print join("\t", @tokenHash{@columnList}), "\n";
		}
		close($reader);
		system("rm $temporaryPrefix.gene.txt");
	}
}

sub writeFastaSequence {
	my ($writer, $sequence) = @_;
	for(my $index = 0; $index < length($sequence); $index += $fastaLineLength) {
		print $writer substr($sequence, $index, $fastaLineLength), "\n";
	}
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
