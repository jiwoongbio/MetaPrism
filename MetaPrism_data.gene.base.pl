#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

(my $codePath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$codePath/data";
system("mkdir -p $dataPath");

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
);
if($help) {
	die <<EOF;

Usage:   perl MetaPrism_data.gene.base.pl [options]

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]

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
my $pid = open2(my $reader, my $writer, 'sort | uniq -c');
forkPrintParentWriter($writer);
if(-r "$dataPath/species_assembly.fasta") {
	system("mkdir -p $dataPath/species_assembly.MetaPrism");
	open(my $reader, "$dataPath/species_assembly.fasta");
	my @speciesIdAssemblyAccessionLineList = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(.*)$/) {
			my ($speciesId, $assemblyAccession, $name) = split(/\|/, $1, 3);
			if(@speciesIdAssemblyAccessionLineList) {
				if($speciesIdAssemblyAccessionLineList[0] != $speciesId || $speciesIdAssemblyAccessionLineList[1] ne $assemblyAccession) {
					if($threads == 1) {
						MetaPrism(@speciesIdAssemblyAccessionLineList);
					} else {
						forkPrintSubroutine(\&MetaPrism, @speciesIdAssemblyAccessionLineList);
					}
					@speciesIdAssemblyAccessionLineList = ();
				}
			}
			($speciesIdAssemblyAccessionLineList[0], $speciesIdAssemblyAccessionLineList[1]) = ($speciesId, $assemblyAccession);
			push(@speciesIdAssemblyAccessionLineList, ">$name");
		} else {
			push(@speciesIdAssemblyAccessionLineList, $line);
		}
	}
	if(@speciesIdAssemblyAccessionLineList) {
		if($threads == 1) {
			MetaPrism(@speciesIdAssemblyAccessionLineList);
		} else {
			forkPrintSubroutine(\&MetaPrism, @speciesIdAssemblyAccessionLineList);
		}
	}
	close($reader);
	forkPrintWait();

	sub MetaPrism {
		my ($speciesId, $assemblyAccession, @lineList) = @_;
		my $file = "$dataPath/species_assembly.MetaPrism/$speciesId.$assemblyAccession.txt";
		{
			(my $logFile = $file) =~ s/\.txt$/.log/;
			open(my $writer, "| perl $codePath/MetaPrism.pl -t $temporaryDirectory -p 1 -m '' - > $file 2> $logFile");
			print $writer "$_\n" foreach(@lineList);
			close($writer);
		}
		if(-s $file) {
			my %geneCountHash = ();
			open(my $reader, $file);
			chomp(my $line = <$reader>);
			my @columnList = split(/\t/, $line, -1);
			while(my $line = <$reader>) {
				chomp($line);
				my %tokenHash = ();
				@tokenHash{@columnList} = split(/\t/, $line, -1);
				$geneCountHash{$tokenHash{'gene'}} += 1;
			}
			close($reader);
			forkPrint("$_\n") foreach(grep {$geneCountHash{$_} == 1} keys %geneCountHash);
		}
	}
}
forkPrintParentWriter();
close($writer);
{
	open(my $writer, "| sort -t '\t' -k2,2nr > $dataPath/single_copy_gene.count.txt");
	while(my $line = <$reader>) {
		chomp($line);
		print $writer "$line\t$1\n" if($line =~ s/^ *([0-9]+) //);
	}
	close($writer);
}
close($reader);
waitpid($pid, 0);
