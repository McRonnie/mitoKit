#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
my $patho_table = "patho_table.txt";
my $MitochondrialDisease = "Mitochondrial-Disease.priority.tab";
my $MutationRNA = "MutationsRNA_lt_MITOMAP_lt_Foswiki.priority.tab";
###
open PATHO, '<', $patho_table;
open CODING, '<', $MitochondrialDisease;
open RNA, '<', $MutationRNA;

my $first_line = <PATHO>;
chomp $first_line;
my @column = split "\t", $first_line;
$column[0] = "conversion";
my $database;
my @order;
while (<PATHO>) {
	chomp;
	my @line = split "\t", $_;
	my $index = 0;
	push @order, $line[0];
	foreach (@column) {
		$database->{$line[0]}->{$column[$index]} = $line[$index];
		$index++;
	}
}
close PATHO;
#print Dumper $database;
while (<CODING>) {
	next if m/del.*ins/;
	next if m/Position/;
	chomp;
	my @line = split "\t", $_;
	my $conversion = $line[3];
	my $transformed = transformConversionStyle($conversion);
	if (exists $database->{$transformed}) {
		$database->{$transformed}->{Priority} = $line[$#line];
	}
	else {
		print STDERR "\e[31mCannot find $conversion:\t$transformed\t$line[$#line]\e[m\n"
	}
}

while (<RNA>) {
	next if m/del.*ins/;
	next if m/Position/;
	chomp;
	my @line = split "\t", $_;
	my $conversion = $line[3];
	my $transformed = transformConversionStyle($conversion);
	if (exists $database->{$transformed}) {
		$database->{$transformed}->{Priority} = $line[$#line];
	}
	else {
		print STDERR "\e[31mCannot find $conversion:\t$transformed\t$line[$#line]\e[m\n"
	}
}
@column = (@column[0..2], "Priority", @column[3..$#column]);
print join "\t", @column;
print "\n";
foreach my $conversion (@order) {
	$database->{$conversion}->{Priority} = 0 unless defined $database->{$conversion}->{Priority};
	print $database->{$conversion}->{conversion};
	foreach my $col (@column[1..$#column]) {
		print "\t", $database->{$conversion}->{$col};
	}
	print "\n";
}



sub transformConversionStyle {
	my $origin = shift @_;
	if ($origin =~ m/(\b[ATCG])(\d+)([ATCG]\b)/) {
		return "$2$3";
	}
	elsif ($origin =~ m/(\d+)del(\d+)/) {
		my $end = $1 + $2 - 1;
		return "${1}-${end}d";
	}
	elsif ($origin =~ m/d/) {
		$origin =~ s/del/d/;
		$origin =~ s/(\d+)//;
		my $first_pos = $1;
		$origin =~ s/(\d+)//;
		my $second_pos = $1;
		if ($origin =~ m/([ATGC]+)/) {
			my $deletion = $1;
			if (length($deletion) == 1) {
				return "${first_pos}d";
			}
			else {
				my $end = $first_pos + length($deletion) - 1;
				return "${first_pos}-${end}d";
			}
		}
	}
	elsif ($origin =~ m/inv/) {
		$origin =~ s/(\d+)//;
		my $first_pos = $1;
		my $inversion = $1 if $origin =~ m/([ATGC]+)/;
		return "$first_pos.$inversion";
	}
	elsif ($origin =~ m/ins/) {
		$origin =~ m/(\d+)ins([ATCG])/;
		my ($pos, $inversion) = ($1, $2);
		return "$pos.$inversion";
	}
	elsif ($origin =~ m/([ATCG]+)(\d+)([ATCG]+)/) {
		my ($left, $pos, $right) = ($1, $2, $3);
		if ($right =~ s/$left//) {
			return "$pos.$right";
		}
		else {
			return "$pos$right"
		}

	}
	else {
		print "Cannot transform $origin\n"
	}
}
