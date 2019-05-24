# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(max);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/MetaPrism_data";
system("mkdir -p $dataPath");

my $database = loadDefaultDatabase();
my $databasePrefix = "$dataPath/$database";

my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	'A=s' => \(my $genomePrefix = ''),
	'B' => \(my $bamInput = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'p=i' => \(my $threads = 1),
	'e=f' => \(my $evalue = 10),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'a=f' => \(my $acceleration = 0.5),
	'C=s' => \@codonList,
	'S=s' => \(my $startCodons = 'GTG,ATG,CTG,TTG,ATA,ATC,ATT'),
	'T=s' => \(my $terminationCodons = 'TAG,TAA,TGA'),
	'l=i' => \(my $minimumTranslationLength = 10),
	'c=f' => \(my $minimumCoverage = 0.8),
	'q=i' => \(my $minimumMappingQuality = 0),
	's=s' => \(my $stranded = ''),
	'P=s' => \(my $contigPrefix = ''),
	'gene2definition=s' => \(my $geneDefinitionFile = ''),
	'protein2gene=s' => \(my $proteinGeneFile = ''),
	'd=s' => \$databasePrefix);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism_gene.pl [options] output.prefix genome.fasta [input.fastq|input.R1.fastq,input.R2.fastq [...]]

Options: -h       display this help message
         -A STR   prepared genome prefix
         -B       input indexed sorted BAM file instead of FASTQ file
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -p INT   number of threads [$threads]
         -e FLOAT maximum e-value to report alignments [$evalue]
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -a FLOAT search acceleration for ublast [$acceleration]
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [$startCodons]
         -T STR   comma-separated termination codons [$terminationCodons]
         -l INT   minimum translation length [$minimumTranslationLength]
         -c FLOAT minimum coverage [$minimumCoverage]
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -s STR   strand specificity, "f" or "r"
         -P STR   contig prefix used for abundance estimation

EOF
}
my ($outputPrefix, $genomeFastaFile, @inputFileList) = @ARGV;
my $genomeNotPrepared = '';
if($genomePrefix eq '') {
	$genomeNotPrepared = 1;
	$genomePrefix = $outputPrefix;
} else {
	chomp(my @fileList = `ls $genomePrefix.bwa_index.*`);
	die "ERROR in $0: '$genomePrefix.bwa_index' is not available.\n" unless(@fileList);
	die "ERROR in $0: '$genomePrefix.region.txt' is not readable.\n" unless(-r "$genomePrefix.region.txt");
}

(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR in $0: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper =~ /^diamond/ || $mapper =~ /^usearch/);
foreach('bwa', 'samtools', $mapperPath) {
	die "ERROR in $0: '$_' is not executable.\n" unless(-x getCommandPath($_));
}
my $databaseFile = '';
$databaseFile = "$databasePrefix.dmnd" if($mapper =~ /^diamond/);
$databaseFile = "$databasePrefix.udb" if($mapper =~ /^usearch/);
foreach($genomeFastaFile, $databaseFile, "$databasePrefix.length.txt") {
	die "ERROR in $0: '$_' is not readable.\n" unless(-r $_);
}
my $outputDirectory = ($outputPrefix =~ /^(.*\/)/) ? $1 : '.';
system("mkdir -p $outputDirectory");
foreach($outputDirectory, $temporaryDirectory) {
	die "ERROR in $0: '$_' is not a writable directory.\n" unless(-d $_ && -w $_);
}

if($databasePrefix eq "$dataPath/$database") {
	if($database =~ /^orthology_uniref/) {
		$geneDefinitionFile = "$dataPath/KEGG_orthology.txt" if($geneDefinitionFile eq '');
	} elsif($database =~ /^ARDB/ || $database =~ /^betalactamases/) {
		$geneDefinitionFile = "$databasePrefix.definition.txt" if($geneDefinitionFile eq '');
		$proteinGeneFile = "$databasePrefix.txt" if($proteinGeneFile eq '');
	}
}

my %geneDefinitionHash = ();
if($geneDefinitionFile ne '') {
	die "ERROR in $0: '$geneDefinitionFile' is not readable.\n" unless(-r $geneDefinitionFile);
	open(my $reader, $geneDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($gene, $definition) = split(/\t/, $line, -1);
		$geneDefinitionHash{$gene} = $definition;
	}
	close($reader);
}
my %proteinGeneHash = ();
my %proteinCutoffHash = ();
if($proteinGeneFile ne '') {
	die "ERROR in $0: '$proteinGeneFile' is not readable.\n" unless(-r $proteinGeneFile);
	open(my $reader, $proteinGeneFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($protein, $gene, $cutoff) = split(/\t/, $line, -1);
		$proteinGeneHash{$protein} = $gene;
		$proteinCutoffHash{$protein} = $cutoff;
	}
	close($reader);
}

if($genomeNotPrepared && $bamInput eq '') { # BWA index
	system("bwa index -p $genomePrefix.bwa_index $genomeFastaFile 1>&2");
}

my @bamFileList = ();
if(@inputFileList) { # Read mapping
	if($bamInput) {
		push(@bamFileList, @inputFileList);
	} else {
		open(my $writer, "| samtools sort -o $outputPrefix.sorted.bam -@ $threads -");
		my $printHeader = 1;
		foreach my $fastqFile (@inputFileList) {
			if($fastqFile =~ /^(.+),(.+)$/) {
				my ($fastqFile1, $fastqFile2) = ($1, $2);
				die "ERROR in $0: '$fastqFile1' is not readable.\n" unless(-r $fastqFile1);
				die "ERROR in $0: '$fastqFile2' is not readable.\n" unless(-r $fastqFile2);
				open(my $reader, "bwa mem -t $threads $genomePrefix.bwa_index $fastqFile1 $fastqFile2 |");
				while(my $line = <$reader>) {
					print $writer $line if($line !~ /^\@/ || $printHeader);
				}
				close($reader);
			} else {
				die "ERROR in $0: '$fastqFile' is not readable.\n" unless(-r $fastqFile);
				open(my $reader, "bwa mem -t $threads $genomePrefix.bwa_index $fastqFile |");
				while(my $line = <$reader>) {
					print $writer $line if($line !~ /^\@/ || $printHeader);
				}
				close($reader);
			}
			$printHeader = '';
		}
		close($writer);
		system("samtools index $outputPrefix.sorted.bam") if(-r "$outputPrefix.sorted.bam");
		if(-r "$outputPrefix.sorted.bam" and -r "$outputPrefix.sorted.bam.bai") {
			chomp(my $readCount = `samtools view -c -F 2304 -q $minimumMappingQuality $outputPrefix.sorted.bam`);
			chomp(my $mappingReadCount = `samtools view -c -F 2308 -q $minimumMappingQuality $outputPrefix.sorted.bam`);
			my $mappingReadRatio = $mappingReadCount / $readCount;
			print "Mapping: $mappingReadCount / $readCount ($mappingReadRatio)\n";
			push(@bamFileList, "$outputPrefix.sorted.bam");
		} else {
			die "ERROR in $0: Read mapping failed.\n";
		}
	}
}

if($genomeNotPrepared) { # ORF translation
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
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. (length($sequence) / 3) - 1);
	}
	my %startCodonHash = map {$_ => 1} split(/,/, $startCodons);
	my %terminationCodonHash = map {$_ => 1} split(/,/, $terminationCodons);
	my $minimumLength = $minimumTranslationLength * 3;
	open(my $reader, $genomeFastaFile);
	open(my $writer, "> $genomePrefix.translation.fasta");
	my ($contig, $sequence) = ('', '');
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			writeTranslationSequences($contig, $sequence) if($contig ne '' && $sequence ne '');
			($contig, $sequence) = ($1, '');
		} else {
			$sequence .= $line;
		}
	}
	writeTranslationSequences($contig, $sequence) if($contig ne '' && $sequence ne '');
	close($reader);
	close($writer);

	sub writeTranslationSequences {
		my ($contig, $sequence) = @_;
		$sequence = uc($sequence);
		my $sequenceLength = length($sequence);
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($sequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $terminationCodonHash{$codon}) {
					writeTranslationSequence($contig, $sequence, $sequenceLength, '+', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
		(my $reverseComplementarySequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($reverseComplementarySequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $terminationCodonHash{$codon}) {
					writeTranslationSequence($contig, $reverseComplementarySequence, $sequenceLength, '-', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
	}

	sub writeTranslationSequence {
		my ($contig, $sequence, $sequenceLength, $strand, $endIndex, @startIndexList) = @_;
		if((my $length = $endIndex - $startIndexList[0]) >= $minimumLength) {
			my ($start, $end) = ('', '');
			($start, $end) = ($startIndexList[0] + 1, $endIndex) if($strand eq '+');
			($start, $end) = (($sequenceLength - $endIndex) + 1, ($sequenceLength - $startIndexList[0])) if($strand eq '-');
			my @startList = map {($_ - $startIndexList[0]) / 3 + 1} @startIndexList;
			print $writer '>', join('|', $contig, $start, $end, $strand, @startList), "\n";
			(my $translationSequence = translate(substr($sequence, $startIndexList[0], $length))) =~ s/\*$//;
			print $writer "$translationSequence\n";
		}
	}
}

my @contigList = ();
if($genomeNotPrepared) { # ORF translation mapping
	system("$mapperPath blastp --query $genomePrefix.translation.fasta --db $databaseFile --out $genomePrefix.blast.txt --outfmt 6 --evalue $evalue --threads $threads --tmpdir $temporaryDirectory --max-target-seqs 0 1>&2") if($mapper =~ /^diamond/);
	system("$mapperPath -ublast $genomePrefix.translation.fasta -db $databaseFile -blast6out $genomePrefix.blast.txt -evalue $evalue -accel $acceleration -threads $threads 1>&2") if($mapper =~ /^usearch/);

	my %proteinSequenceLengthHash = ();
	{
		open(my $reader, "$databasePrefix.length.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my ($protein, $sequenceLength) = split(/\t/, $line, -1);
			$proteinSequenceLengthHash{$protein} = $sequenceLength;
		}
		close($reader);
	}
	{
		open(my $reader, "$genomePrefix.blast.txt");
		open(my $writer, "> $genomePrefix.blast.filtered.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)} = split(/\t/, $line, -1);
			next if(defined($_ = $proteinCutoffHash{$tokenHash{'sseqid'}}) && $tokenHash{'pident'} < $_);
			$tokenHash{'coverage'} = (($tokenHash{'pident'} / 100) * $tokenHash{'length'}) / $proteinSequenceLengthHash{$tokenHash{'sseqid'}};
			print $writer join("\t", $line, $tokenHash{'coverage'}), "\n" if($tokenHash{'coverage'} >= $minimumCoverage && $tokenHash{'qstart'} <= $tokenHash{'qend'} && $tokenHash{'sstart'} <= $tokenHash{'send'});
		}
		close($reader);
		close($writer);
	}
	{
		my @columnList = ('contig', 'start', 'end', 'strand', 'gene');
		push(@columnList, 'definition') if(%geneDefinitionHash);
		{
			open(my $writer, "> $genomePrefix.region.txt");
			print $writer join("\t", @columnList), "\n";
			close($writer);
		}
		my %contigHash = ();
		open(my $reader, "sort -t '\t' -k1,1 -k11,11g -k12,12gr -k13,13gr $genomePrefix.blast.filtered.txt |");
		open(my $writer, "| sort -t '\t' -k1,1 -k2,2n -k3,3n -k4,4r -k5,5 | uniq >> $genomePrefix.region.txt");
		my %topTokenHash = ();
		$topTokenHash{'qseqid'} = '';
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore), 'coverage'} = split(/\t/, $line, -1);
			%topTokenHash = %tokenHash if($tokenHash{'qseqid'} ne $topTokenHash{'qseqid'});
			if(scalar(grep {$tokenHash{$_} != $topTokenHash{$_}} 'evalue', 'bitscore', 'coverage') == 0) {
				(@tokenHash{'contig', 'start', 'end', 'strand'}, my @startList) = split(/\|/, $tokenHash{'qseqid'});
				if($tokenHash{'qstart'} > 1) {
					my $startIndex = max(grep {$_ <= $tokenHash{'qstart'}} @startList) - 1;
					$tokenHash{'start'} += $startIndex * 3 if($tokenHash{'strand'} eq '+');
					$tokenHash{'end'}   -= $startIndex * 3 if($tokenHash{'strand'} eq '-');
				}
				$tokenHash{'gene'} = getGene($tokenHash{'sseqid'});
				next if($tokenHash{'gene'} eq 'off-target');
				$tokenHash{'definition'} = $geneDefinitionHash{$tokenHash{'gene'}};
				$tokenHash{'definition'} = '' unless(defined($tokenHash{'definition'}));
				print $writer join("\t", @tokenHash{@columnList}), "\n";
				$contigHash{$tokenHash{'contig'}} = 1;
			}
		}
		close($reader);
		close($writer);
		@contigList = sort keys %contigHash;
	}
}

if(@bamFileList) { # Abundance estimation
	my %contigSequenceLengthHash = ();
	{
		open(my $reader, $genomeFastaFile);
		my $contig = '';
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ /^>(\S*)/) {
				$contig = $1;
			} else {
				$contigSequenceLengthHash{$contig} += length($line);
			}
		}
		close($reader);
	}
	my %contigMeanDepthHash = ();
	my $genomeMeanDepth = 0;
	{
		foreach my $bamFile (@bamFileList) {
			foreach my $contig (@contigList) {
				open(my $reader, "samtools view -F 2308 -q $minimumMappingQuality $bamFile $contig |");
				while(my $line = <$reader>) {
					chomp($line);
					my %tokenHash = ();
					@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line, -1);
					my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
					if(@positionList = grep {defined} @positionList) {
						$contigMeanDepthHash{$contig} += scalar(@positionList);
						$genomeMeanDepth += scalar(@positionList);
					}
				}
				close($reader);
			}
		}
		my $genomeSequenceLength = 0;
		foreach my $contig (keys %contigMeanDepthHash) {
			my $contigSequenceLength = $contigSequenceLengthHash{$contig};
			$contigMeanDepthHash{$contig} /= $contigSequenceLength;
			$genomeSequenceLength += $contigSequenceLength;
		}
		$genomeMeanDepth /= $genomeSequenceLength;
		print "Genome mean depth: $genomeMeanDepth\n";
	}
	my @valueColumnList = ('readCount', 'RPKM', 'meanDepth', 'meanDepth/contig', 'meanDepth/genome');
	my %geneTokenHashHash = ();
	{
		my %readCountHash = ();
		my @tokenHashList = ();
		open(my $reader, "$genomePrefix.region.txt");
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			next unless($tokenHash{'contig'} =~ /^$contigPrefix/);
			foreach my $bamFile (@bamFileList) {
				my ($baseCount, @readList) = getBaseCountReadList($bamFile, @tokenHash{'contig', 'start', 'end', 'strand'});
				$tokenHash{'baseCount'} += $baseCount;
				$tokenHash{'readCount'} += scalar(@readList);
				$readCountHash{$_} += 1 foreach(@readList);
			}
			$tokenHash{'length'} = $tokenHash{'end'} - $tokenHash{'start'} + 1;
			$tokenHash{'RPK'} = $tokenHash{'readCount'} / ($tokenHash{'length'} / 1000);
			if(($tokenHash{'meanDepth'} = $tokenHash{'baseCount'} / $tokenHash{'length'}) == 0) {
				$tokenHash{'meanDepth/contig'} = 0;
				$tokenHash{'meanDepth/genome'} = 0;
			} else {
				$tokenHash{'meanDepth/contig'} = $tokenHash{'meanDepth'} / $contigMeanDepthHash{$tokenHash{'contig'}};
				$tokenHash{'meanDepth/genome'} = $tokenHash{'meanDepth'} / $genomeMeanDepth;
			}
			push(@tokenHashList, \%tokenHash);
		}
		close($reader);
		my $totalReadCount = scalar(keys %readCountHash);
		open(my $writer, "> $outputPrefix.region.abundance.txt");
		print $writer join("\t", @columnList, @valueColumnList), "\n";
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			$tokenHash{'RPKM'} = $tokenHash{'RPK'} / ($totalReadCount / 1000000);
			print $writer join("\t", @tokenHash{@columnList, @valueColumnList}), "\n";
			$geneTokenHashHash{$tokenHash{'gene'}}->{$_} += $tokenHash{$_} foreach(@valueColumnList);
		}
		close($writer);
	}
	{
		open(my $writer, "> $outputPrefix.abundance.txt");
		if(%geneDefinitionHash) {
			print $writer join("\t", 'gene', 'definition', @valueColumnList), "\n";
			foreach my $gene (sort keys %geneTokenHashHash) {
				my $tokenHash = $geneTokenHashHash{$gene};
				my $definition = $geneDefinitionHash{$gene};
				$definition = '' unless(defined($definition));
				print $writer join("\t", $gene, $definition, @$tokenHash{@valueColumnList}), "\n";
			}
		} else {
			print $writer join("\t", 'gene', @valueColumnList), "\n";
			foreach my $gene (sort keys %geneTokenHashHash) {
				my $tokenHash = $geneTokenHashHash{$gene};
				print $writer join("\t", $gene, @$tokenHash{@valueColumnList}), "\n";
			}
		}
		close($writer);
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

sub getGene {
	my ($protein) = @_;
	return $_ if(defined($_ = $proteinGeneHash{$protein}));
	return $1 if($protein =~ /^(K[0-9]{5})_/);
	return $protein;
}

sub getBaseCountReadList {
	my ($bamFile, $contig, $start, $end, $strand) = @_;
	my $baseCount = 0;
	my %readCountHash = ();
	open(my $reader, "samtools view -F 2308 -q $minimumMappingQuality $bamFile $contig:$start-$end |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line, -1);
		next if($stranded ne '' && getReadStrand($tokenHash{'flag'}) ne $strand);
		my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
		if(@positionList = grep {$start <= $_ && $_ <= $end} grep {defined} @positionList) {
			$baseCount += scalar(@positionList);
			$readCountHash{$tokenHash{'qname'}} += 1;
		}
	}
	close($reader);
	return ($baseCount, keys %readCountHash);
}

sub getReadStrand {
	my ($flag) = @_;
	if($stranded eq 'f') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (99, 147));
			return '-' if(grep {$_ == $flag} (83, 163));
		} else {
			return '+' if($flag == 0);
			return '-' if($flag == 16);
		}
	}
	if($stranded eq 'r') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (83, 163));
			return '-' if(grep {$_ == $flag} (99, 147));
		} else {
			return '+' if($flag == 16);
			return '-' if($flag == 0);
		}
	}
	return '';
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		}
		if($operation eq 'I') {
			$index += $length;
		}
		if($operation eq 'D') {
			$position += $length;
		}
		if($operation eq 'N') {
			$position += $length;
		}
		if($operation eq 'S') {
			$index += $length;
		}
	}
	return @positionList;
}
