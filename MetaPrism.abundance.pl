#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use IPC::Open2;
use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

(my $code_path = abs_path($0)) =~ s/\/[^\/]*$//;
my $data_path = "$code_path/data";
checked_system('mkdir', '-p', $data_path);

my @sam_mandatory_column_list = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');

my $default_base_abundance_genes = 'K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965';
my $default_number_of_base_abundance_genes = 10;

chomp(my $hostname = `hostname`);
$? == 0 or die "hostname failed (status $?)\n";
my $temporary_directory = $ENV{'TMPDIR'};
$temporary_directory = '/tmp' unless($temporary_directory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporary_directory,
	'p=i' => \(my $threads = 1),
	'number_per_thread=i' => \(my $number_per_thread = 1000),
	'q=i' => \(my $minimum_mapping_quality = 0),
	'f=i' => \(my $include_flag = 0),
	'F=i' => \(my $exclude_flag = 0),
	'S=s' => \(my $stranded = ''),
	'B=f' => \(my $base_abundance = ''),
	'b=s' => \(my $base_abundance_genes),
	'base_abundance_output_file=s' => \(my $base_abundance_output_file = ''),
	'use_merged_bam_file' => \(my $use_merged_bam_file = ''),
	'use_samtools_depth' => \(my $use_samtools_depth = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl MetaPrism.abundance.pl [options] MetaPrism.txt genome.sorted.bam [...] > MetaPrism.abundance.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -q INT   minimum mapping quality [$minimum_mapping_quality]
         -f INT   include flag [$include_flag]
         -F INT   exclude flag [$exclude_flag]
         -S STR   stranded, "f" or "r"
         -B FLOAT base abundance
         -b STR   base abundance genes or number of base abundance genes [$default_base_abundance_genes or $default_number_of_base_abundance_genes if $data_path/single_copy_gene.count.txt is available]
         --base_abundance_output_file
         --use_merged_bam_file use a temporary merged bam file
         --use_samtools_depth use samtools depth instead of parsing cigar strings

EOF
}

$threads >= 1 or die "Number of threads must be at least 1\n";
($stranded eq '' || $stranded eq 'f' || $stranded eq 'r' || $stranded eq 'forward' || $stranded eq 'reverse')
	or die "Invalid -S value '$stranded': expected f, r, forward, or reverse\n";

my ($metaPrismFile, @bamFileList) = @ARGV;
@bamFileList or die "At least one BAM file is required\n";
foreach my $bamFile (@bamFileList) {
	(-f $bamFile && -r $bamFile) or die "Can't read BAM file '$bamFile'\n";
}

my $temporary_prefix = "$temporary_directory/$hostname.$$";
system("rm -f $temporary_prefix.*") == 0 or die "Can't remove '$temporary_prefix.*'\n";
{
	my %pid_hash = ();
	my $writer;
	my $parent_writer;

	sub fork_print_parent_writer {
		($parent_writer) = @_;
	}

	sub fork_print_subroutine {
		my ($subroutine, @arguments) = @_;
		my $pid = fork();
		die "Can't fork: $!" unless(defined($pid));
		if($pid) {
			$pid_hash{$pid} = 1;
		} else {
			open($writer, '>', "$temporary_prefix.$$") or die "Can't open '$temporary_prefix.$$': $!";
			$subroutine->(@arguments);
			close($writer) or die "Can't close '$temporary_prefix.$$': $!";
			exit(0);
		}
		fork_print_wait($threads);
	}

	sub fork_print_wait {
		my ($number) = (@_, 1);
		while(scalar(keys %pid_hash) >= $number) {
			my $pid = wait();
			die "wait failed: $!\n" if($pid == -1);
			die "Child process '$pid' failed\n" if($? != 0);
			if($pid_hash{$pid}) {
				open(my $reader, '<', "$temporary_prefix.$pid") or die "Can't open '$temporary_prefix.$pid': $!";
				if(defined($parent_writer)) {
					print $parent_writer $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader) or die "Can't close '$temporary_prefix.$pid': $!";
				unlink("$temporary_prefix.$pid") or die "Can't remove '$temporary_prefix.$pid': $!";
				delete $pid_hash{$pid};
			}
		}
	}

	sub fork_print {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parent_writer)) {
			print $parent_writer @_;
		} else {
			print @_;
		}
	}
}

if($use_merged_bam_file || $use_samtools_depth) {
	if(scalar(@bamFileList) == 1) {
		checked_system('ln', '-sfr', $bamFileList[0], "$temporary_prefix.bam");
	} else {
		checked_system('samtools', 'merge', '--threads', $threads, "$temporary_prefix.bam", @bamFileList);
	}
	if($use_samtools_depth) {
		if($minimum_mapping_quality == 0 && $include_flag == 0 && $exclude_flag == 0) {
			checked_system('ln', '-sfr', "$temporary_prefix.bam", "$temporary_prefix.filtered.bam");
		} else {
			checked_system("samtools view -b -q $minimum_mapping_quality -f $include_flag -F $exclude_flag $temporary_prefix.bam > $temporary_prefix.filtered.bam");
		}
		if($stranded) {
			open(my $reader, '-|', 'samtools', 'view', '-h', "$temporary_prefix.filtered.bam") or die "Can't run samtools view: $!\n";
			open(my $writerF, '|-', 'samtools', 'view', '-b', '-o', "$temporary_prefix.filtered.f.bam", '-') or die "Can't run samtools view (f): $!\n";
			open(my $writerR, '|-', 'samtools', 'view', '-b', '-o', "$temporary_prefix.filtered.r.bam", '-') or die "Can't run samtools view (r): $!\n";
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ /^\@/) {
					print $writerF $line, "\n";
					print $writerR $line, "\n";
					next;
				}
				my %token_hash = ();
				(@token_hash{@sam_mandatory_column_list}, my @tagTypeValueList) = split(/\t/, $line, -1);
				$token_hash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
				print $writerF $line, "\n" if(grep {($token_hash{'flag'} & 253) == $_} 97, 145, 0);
				print $writerR $line, "\n" if(grep {($token_hash{'flag'} & 253) == $_} 81, 161, 16);
			}
			close($reader) or die "samtools view failed (status $?)\n";
			close($writerF) or die "samtools view (f) failed (status $?)\n";
			close($writerR) or die "samtools view (r) failed (status $?)\n";
			checked_system('samtools', 'index', "$temporary_prefix.filtered.f.bam");
			checked_system('samtools', 'index', "$temporary_prefix.filtered.r.bam");
		} else {
			checked_system('samtools', 'index', "$temporary_prefix.filtered.bam");
		}
	} else {
		checked_system('samtools', 'index', "$temporary_prefix.bam");
		@bamFileList = ("$temporary_prefix.bam");
	}
}
my $use_sam_module = !$use_samtools_depth && eval {require Bio::DB::Sam; 1;};
unless($use_samtools_depth) {
	foreach my $bamFile (@bamFileList) {
		next if(-r "$bamFile.bai");
		(my $baiFile = $bamFile) =~ s/\.bam$/.bai/;
		if(-r $baiFile) {
			checked_system('ln', '-sfr', $baiFile, "$bamFile.bai");
		} else {
			checked_system('samtools', 'index', $bamFile);
		}
	}
}
my $pid = open2(my $reader, my $writer, 'sort -n');
fork_print_parent_writer($writer);
my @column_list = ();
my $gene_index;
{
	open(my $reader, '<', $metaPrismFile) or die "Can't open '$metaPrismFile': $!";
	chomp(my $line = <$reader>);
	$line =~ s/\r$//;
	@column_list = split(/\t/, $line, -1);
	($gene_index) = grep {$column_list[$_] eq 'gene'} 0 .. $#column_list;
	if($base_abundance ne '') {
		$base_abundance_genes = '';
	} elsif(!defined($base_abundance_genes) && defined($gene_index)) {
		$base_abundance_genes = $default_base_abundance_genes;
		$base_abundance_genes = $default_number_of_base_abundance_genes if(-s "$data_path/single_copy_gene.count.txt" and -r "$data_path/single_copy_gene.count.txt");
	}
	if(defined($base_abundance_genes) && $base_abundance_genes ne '') {
		die "Missing 'gene' column required by -b\n" if(!defined($gene_index));
		if(-r $base_abundance_genes) {
			my @base_abundance_gene_list = ();
			open(my $reader, '<', $base_abundance_genes) or die "Can't open '$base_abundance_genes': $!";
			while(my $line = <$reader>) {
				chomp($line);
				push(@base_abundance_gene_list, $line);
			}
			close($reader) or die "Can't close '$base_abundance_genes': $!";
			$base_abundance_genes = join(',', @base_abundance_gene_list);
		} elsif($base_abundance_genes =~ /^[0-9]+$/) {
			my @base_abundance_gene_list = ();
			open(my $reader, '<', "$data_path/single_copy_gene.count.txt") or die "Can't open '$data_path/single_copy_gene.count.txt': $!";
			while(scalar(@base_abundance_gene_list) < $base_abundance_genes) {
				chomp(my $line = <$reader>);
				my ($base_abundance_gene) = split(/\t/, $line, -1);
				push(@base_abundance_gene_list, $base_abundance_gene);
			}
			close($reader) or die "Can't close '$data_path/single_copy_gene.count.txt': $!";
			$base_abundance_genes = join(',', @base_abundance_gene_list);
		}
	}
	my $order = 0;
	my @token_list_list = ();
	while(my $line = <$reader>) {
		chomp($line);
		$line =~ s/\r$//;
		my @token_list = split(/\t/, $line, -1);
		push(@token_list_list, [$order += 1, @token_list]);
		if(scalar(@token_list_list) >= $number_per_thread) {
			if($threads == 1) {
				printAbundance(@token_list_list);
			} else {
				fork_print_subroutine(\&printAbundance, @token_list_list);
			}
			@token_list_list = ();
		}
	}
	if(@token_list_list) {
		if($threads == 1) {
			printAbundance(@token_list_list);
		} else {
			fork_print_subroutine(\&printAbundance, @token_list_list);
		}
	}
	close($reader) or die "Can't close '$metaPrismFile': $!";
	fork_print_wait();
}
fork_print_parent_writer();
close($writer) or die "Can't close sort input: $!\n";
{
	my @token_list_abundance_list = ();
	my %gene_abundance_hash = ();
	while(my $line = <$reader>) {
		chomp($line);
		my ($order, @token_list) = split(/\t/, $line, -1);
		my $abundance = pop(@token_list);
		push(@token_list_abundance_list, [\@token_list, $abundance]);
		$gene_abundance_hash{$token_list[$gene_index]} += $abundance if(defined($gene_index));
	}
	if(defined($base_abundance_genes) && $base_abundance_genes ne '') {
		my @abundanceList = map {defined($_) ? $_ : 0} @gene_abundance_hash{split(/,/, $base_abundance_genes)};
		$base_abundance = median(@abundanceList);
		if($base_abundance_output_file ne '') {
			open(my $writer, '>', $base_abundance_output_file) or die "Can't open '$base_abundance_output_file': $!";
			print $writer $base_abundance, "\n";
			close($writer) or die "Can't close '$base_abundance_output_file': $!";
		}
	}
	print join("\t", @column_list, 'abundance'), "\n";
	if($base_abundance eq '') {
		foreach(@token_list_abundance_list) {
			my ($token_list, $abundance) = @$_;
			print join("\t", @$token_list, $abundance), "\n";
		}
	} elsif($base_abundance > 0) {
		foreach(@token_list_abundance_list) {
			my ($token_list, $abundance) = @$_;
			print join("\t", @$token_list, $abundance / $base_abundance), "\n";
		}
	} else {
		foreach(@token_list_abundance_list) {
			my ($token_list, $abundance) = @$_;
			print join("\t", @$token_list, "$abundance/$base_abundance"), "\n";
		}
	}
}
close($reader) or die "Can't close sort output: $!\n";
waitpid($pid, 0) == $pid or die "waitpid for sort failed: $!\n";
$? == 0 or die "sort failed (status $?)\n";
if($use_merged_bam_file || $use_samtools_depth) {
	checked_system("rm $temporary_prefix.bam");
	if($use_samtools_depth) {
		checked_system("rm $temporary_prefix.filtered.bam");
		if($stranded) {
			checked_system("rm $temporary_prefix.filtered.f.bam");
			checked_system("rm $temporary_prefix.filtered.f.bam.bai");
			checked_system("rm $temporary_prefix.filtered.r.bam");
			checked_system("rm $temporary_prefix.filtered.r.bam.bai");
		} else {
			checked_system("rm $temporary_prefix.filtered.bam.bai");
		}
	} else {
		checked_system("rm $temporary_prefix.bam.bai");
	}
}

sub printAbundance {
	my @samList = map {Bio::DB::Sam->new(-bam => $_)} @bamFileList if($use_sam_module);
	foreach(@_) {
		my ($order, @token_list) = @$_;
		my %token_hash = ();
		@token_hash{@column_list} = @token_list;
		my ($chromosome, $start, $end, $strand);
		($chromosome, $start, $end, $strand) = @token_hash{'chromosome', 'start', 'end', 'strand'} if(all {defined} @token_hash{'chromosome', 'start', 'end', 'strand'});
		($chromosome, $start, $end, $strand) = @token_hash{'Contig id', 'Start', 'Stop', 'Strand'} if(all {defined} @token_hash{'Contig id', 'Start', 'Stop', 'Strand'});
		my $depth = 0;
		if($use_samtools_depth) {
			if($stranded eq '') {
				$depth = get_bam_file_region_depth("$temporary_prefix.filtered.bam", $chromosome, $start, $end);
			} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
				$depth = get_bam_file_region_depth("$temporary_prefix.filtered.f.bam", $chromosome, $start, $end);
			} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
				$depth = get_bam_file_region_depth("$temporary_prefix.filtered.r.bam", $chromosome, $start, $end);
			}
		} elsif($use_sam_module) {
			foreach my $sam (@samList) {
				foreach my $alignment ($sam->get_features_by_location(-seq_id => $chromosome, -start => $start, -end => $end)) {
					next if($alignment->qual < $minimum_mapping_quality);
					next if(($alignment->flag & int($include_flag)) != $include_flag);
					next if($alignment->flag & int($exclude_flag));
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($alignment->flag & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($alignment->flag & 253) == $_} 81, 161, 16);
					}
					my @position_list = get_position_list($alignment->start, $alignment->cigar_str);
					@position_list = grep {$_ ne ''} @position_list;
					@position_list = grep {$start <= $_ && $_ <= $end} @position_list;
					$depth += scalar(@position_list);
				}
			}
		} else {
			foreach my $bamFile (@bamFileList) {
				my $command;
				if($stranded eq '') {
					$command = "samtools view -h -q $minimum_mapping_quality -f $include_flag -F $exclude_flag $bamFile $chromosome:$start-$end | samtools depth -d 0 -";
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
					$command = "samtools view -h -q $minimum_mapping_quality -f $include_flag -F $exclude_flag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@token_list = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$token_list[1] & 253) == \$_} 97, 145, 0);' | samtools depth -d 0 -";
				} elsif((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
					$command = "samtools view -h -q $minimum_mapping_quality -f $include_flag -F $exclude_flag $bamFile $chromosome:$start-$end | perl -Mstrict -Mwarnings -wne 'chomp(my \$line = \$_); my \@token_list = split(/\\t/, \$line, -1); print \"\$line\\n\" if(\$line =~ /^\@/ || grep {(\$token_list[1] & 253) == \$_} 81, 161, 16);' | samtools depth -d 0 -";
				}
				open(my $reader, '-|', 'bash', '-o', 'pipefail', '-c', $command) or die "Can't run depth pipeline: $!\n";
				while(my $line = <$reader>) {
					chomp($line);
					my ($chromosome, $position, $depth1) = split(/\t/, $line, -1);
					$depth += $depth1 if($start <= $position && $position <= $end);
				}
				close($reader) or die "Depth pipeline failed (status $? ): $command\n";
			}
		}
		my $length = $end - $start + 1;
		my $meanDepth = $depth / $length;
		my $abundance = $meanDepth;
		fork_print(join("\t", $order, @token_list, $abundance), "\n");
	}
}

sub get_bam_file_region_depth {
	my ($bamFile, $chromosome, $start, $end) = @_;
	my $region_depth = 0;
	open(my $reader, '-|', 'samtools', 'depth', '-a', '-a', '-r', "$chromosome:$start-$end", $bamFile) or die "Can't run samtools depth: $!\n";
	while(my $line = <$reader>) {
		chomp($line);
		my ($chromosome, $position, $position_depth) = split(/\t/, $line, -1);
		$region_depth += $position_depth;
	}
	close($reader) or die "samtools depth failed (status $?)\n";
	return $region_depth;
}

sub checked_system {
	system(@_) == 0 or die "Command failed (status $? ): @_\n";
}

sub get_position_list {
	my ($position, $cigar) = @_;
	my @position_list = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M' || $operation eq '=' || $operation eq 'X') {
			@position_list[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I' || $operation eq 'S') {
			@position_list[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D' || $operation eq 'N') {
			$position += $length;
		}
	}
	return @position_list;
}

sub median {
	my @token_list = sort {$a <=> $b} @_;
	return undef unless(@token_list);
	my $number = scalar(@token_list);
	if($number % 2 == 0) {
		return ($token_list[$number / 2 - 1] + $token_list[$number / 2]) / 2;
	} else {
		return $token_list[($number - 1) / 2];
	}
}
