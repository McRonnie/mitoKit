#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/max min/;
use Text::NSP::Measures::2D::Fisher::twotailed;


##Usage



#Specify final file from MToolBox and child tab file from mitoKit
my $path;
my $debug = 0;
GetOptions(
		   'final|f=s' => \ my $finalPath,
		   'childp|c=s' => \ $path->{childp},
		   'mother|m=s' => \ $path->{mother},
		   'debug|d' => \ $debug,
		  );
my $parsedData;
for my $childOrMother (qw/childp mother/) {
	open my $fl, '<', $path->{$childOrMother};
	my @header;
	while (<$fl>) {
		chomp;
		if (/^#/) {
			s/#//;
			@header = split "\t", $_;
			next;
		}
		my @line = split "\t", $_;
		my $pos = $line[1];
		for my $col (0..$#header) {
			$parsedData->{$childOrMother}->{$pos}->{$header[$col]} = $line[$col] if defined $line[$col];
		}
	}
	close $fl;
}

#print Dumper $parsedData if $debug;

my $tmp;
map {$tmp->{$_} = $_} qw/A T G C/;
##Precessing <final> file from MToolBox
open my $finalFL, '<', $finalPath;
my @finalHeader;
my ($childpCol, $motherCol, $posCol) = (0, 0, 0);
print join("\t", "#Pos", "Sample", "A(child|mother)", "T(child|mother)", "G(child|mother)", "C(child|mother)", "Major/Minor_Base_Consistency", "Fisher_Exact_P_Value"), "\n";
while (<$finalFL>) {
	chomp;
	if (/^Child/) {
		@finalHeader = split "\t", $_;
		$childpCol++ while !($finalHeader[$childpCol] =~ m/Child/);		##Catch column number of child and mother from <final> file
		$motherCol++ while !($finalHeader[$motherCol] =~ m/mother/);
		$posCol++ while !($finalHeader[$posCol] =~ m/Variant Allel/);
#		print join("\n", $childpCol, $motherCol, $posCol), "\n";
		next;
	}
	my @line = split "\t", $_;
	my ($childpBaseCount, $motherBaseCount, $pos) = @line[$childpCol, $motherCol, $posCol];
	my $finalDataset;
	##Catch information in <final> file
	my @modelBaseCount = qw/Pos Ref Alt Depth MeanQ A C G T/; ##For easy use
	$childpBaseCount =~ m/(\d+) ([ATGC]) ([ATGC]) (\d+) ([\d.]+) \((\d+), (\d+), (\d+), (\d+)\)/;
	for my $i (0..$#modelBaseCount) {
		$finalDataset->{childp}->{$modelBaseCount[$i]} = ($1, $2, $3, $4, $5, $6, $7, $8, $9)[$i];
	}
	$motherBaseCount =~ m/(\d+) ([ATGC]) ([ATGC]) (\d+) ([\d.]+) \((\d+), (\d+), (\d+), (\d+)\)/;
	for my $i (0..$#modelBaseCount) {
		$finalDataset->{mother}->{$modelBaseCount[$i]} = ($1, $2, $3, $4, $5, $6, $7, $8, $9)[$i];
	}
	##Parameters used for comparison
	my $mafProfiling = 'Unmatched';
	my $chisqP = {childp => 'NA',
				  mother => 'NA'};

	foreach my $childOrMother (qw/childp mother/) {
		print STDERR join("\n", map {$parsedData->{$childOrMother}->{$pos}->{$_}} qw/A T G C/), "\n" if $debug;
		my ($parsedMajor, $parsedMinor) =
			(
			 sort
			 {
				 $parsedData->{$childOrMother}->{$pos}->{$b}
					 <=>
					 $parsedData->{$childOrMother}->{$pos}->{$a};
			 }
			 qw/A T G C/
			);
		my ($finalMajor, $finalMinor) =
			(
			 sort
			 {
				 $finalDataset->{$childOrMother}->{$b}
					 <=>
					 $finalDataset->{$childOrMother}->{$a};
			 }
			 qw/A T G C/
			)[0,1];
		if ($parsedMajor eq $finalMajor and $parsedMinor eq $finalMinor) {
			$mafProfiling = 'Matched';
		}
		$chisqP->{$childOrMother} = fisher(
										   $parsedData->{$childOrMother}->{$pos}->{$parsedMajor},
										   $parsedData->{$childOrMother}->{$pos}->{$parsedMinor},
										   $finalDataset->{$childOrMother}->{$finalMajor},
										   $finalDataset->{$childOrMother}->{$finalMinor},
										  );
		print STDERR Dumper $finalDataset if $debug;
		print join("\t",
				   (
					$pos,
					$childOrMother eq 'mother' ? "mother" : "child",
					map({
						 $parsedData->{$childOrMother}->{$pos}->{$_}."|".
						 $finalDataset->{$childOrMother}->{$_}
						} qw/A T G C/),
					$mafProfiling,
					$chisqP->{$childOrMother},
				   )
				  ),"\n";
	}
}






sub fisher {
	my ($n11, $n12, $n21, $n22) = @_;
	my $n1p = $n11 + $n12;
	my $np1 = $n11 + $n21;
	my $npp = $n1p + $n21 + $n22;
	calculateStatistic(n11 => $n11,
					   n1p => $n1p,
					   np1 => $np1,
					   npp => $npp
					  );
}
