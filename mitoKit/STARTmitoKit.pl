#!/usr/bin/env perl
use strict;
use warnings;
use Mito;
use File::Basename;
use File::Path qw/make_path/;
use Getopt::Long;
use File::Copy;
our $VERSION = 0.1;

##Flags for jumping
my $debug = 0;
my $anno = 0;
my $pileup = 0;
GetOptions(
		   "debug" => \ $debug,
		   "anno" => \ $anno,
		   "pileup" => \ $pileup,
		  );







##Start
open SAMPLE, '<', "SAMPLE"
    or die "\e[31mCannot open <SAMPLE> file\e[m";
open CONFIG, '<', "CONFIG"
    or die "\e[31mCannot open <CONFIG> file\e[m";
our ($workdir, $childFastq1, $childFastq2, $childFull1, $childFull2, $motherFastq1, $motherFastq2, $motherFull1, $motherFull2, $childBase, $motherBase);
##our $config_parameters;

our ($threads, $threads_0, $chrM_ref, $node_name, $picard_path);

##Read work directory and file names from <SAMPLE> file
while (<SAMPLE>) {
    chomp;
    if (/Workdir/) {
	$workdir = $1 if m/Workdir:*\s*(\S+)$/;
	$workdir .= "/" unless $workdir =~ m|/$|;
    }
    elsif (/Child/) {
	m/Child:*\s*(\S+)\s*(\S+)/;
	($childFastq1, $childFastq2) = ($1, $2);
	($childFull1, $childFull2) = map {$workdir . "0_Raw/" . $_} ($childFastq1, $childFastq2);
	$childFastq1 =~ m/([\._]R[12].*fastq.*$)/;
	my $suffix = $1;
	$childBase = $childFastq1 =~ s/$suffix//r;
    }
    elsif (/Mother/) {
	m/Mother:*\s*(\S+)\s*(\S+)/;
	($motherFastq1, $motherFastq2) = ($1, $2);
	($motherFull1, $motherFull2) = map {$workdir . "0_Raw/" . $_} ($motherFastq1, $motherFastq2);
	$motherFastq1 =~ m/([\._]R[12].*fastq.*$)/;
	my $suffix = $1;
	$motherBase = $motherFastq1 =~ s/$suffix//r;
    }
}
close SAMPLE;


##Read reference file path from <CONFIG> file
while (<CONFIG>) {
    next if m/#/;
    next if m/^\s*$/;
    chomp;
    if (m/chrM_ref/) {
	$chrM_ref = ReadCONFIG($_);
    }
    elsif (m/threads/) {
	$threads = ReadCONFIG($_);
	$threads_0 = $threads - 1;
    }
    elsif (m/node/) {
	$node_name = ReadCONFIG($_);
    }
    elsif (m/picard/) {
	$picard_path = ReadCONFIG($_);
    }
#   else {
#	m/(\S+)\s*=\s*(\S+)/;
#	$config_parameters->{$1} = $2;
#   }
}
##Set parameters used in Pileup.pm
#our ($chrM_bed, $chrM_ref, $table_path, $min_depth, $hetero_threshold) =
#    map {$config_parameters->{$_}} qw/chrM_bed chrM_ref table_path min_depth hetero_threshold/;

goto ANNO if $anno;
goto PILEUP if $pileup;

###START###
chdir $workdir;
##write standard error info into $workdir/log file
#open STDERR, '>', "log";

##Generate sub.sh

=DEPRECATED
open SUB, '>', "sub.sh";
print SUB
    "#!/usr/bin/sh
#SBATCH -p $node_name
#SBATCH -J $childBase
#SBATCH -N 1
#SBATCH -n $threads
";
=cut

unless (-e "1_Mapping") {
    make_path("1_Mapping") or die "\e[31mCannot make path ${workdir}1_Mapping\e[m";
    print STDERR "\e[32mMake directory ${workdir}1_Mapping\e[m\n";
}


unless (-e "2_Pileup") {
    make_path("2_Pileup") or die "\e[31mCannot make path ${workdir}2_Pileup\e[m";
    print STDERR "\e[32mMake directory ${workdir}2_Pileup\e[m\n";
}


unless (-e "3_Annotation") {
    make_path("3_Annotation") or die "\e[31mCannot make path ${workdir}3_Annotation\e[m";
    print STDERR "\e[32mMake directory ${workdir}3_Annotation\e[m\n";
}



##Mapping to reference using bwa mem
###Child
chdir($workdir . "1_Mapping") or die "\e[31mCannot enter ${workdir}1_Mapping\e[m";
print STDERR "\e[32mChange current directory to ${workdir}1_Mapping\e[m\n";

=DEPRECATED
print SUB "cd ${workdir}1_Mapping
";
print SUB "bwa mem -t $threads $chrM_ref $childFull1 $childFull2 -R \'\@RG\\tID:$childBase\\tSM:$childBase\\tPL:ILLUMINA\' | samtools view -@ $threads_0 -O bam -o $childBase.bam &&
";
print SUB "samtools fixmate -@ $threads_0 -O bam $childBase.bam $childBase.fixmate.bam &&
";
print SUB "samtools sort -@ $threads_0 -O bam -o $childBase.sorted.bam $childBase.fixmate.bam &&
";
###Mother
print SUB "bwa mem -t $threads $chrM_ref $motherFull1 $motherFull2 -R \'\@RG\\tID:$motherBase\\tSM:$motherBase\\tPL:ILLUMINA\' | samtools view -@ $threads_0 -O bam -o $motherBase.bam &&
";
print SUB "samtools fixmate -@ $threads_0 -O bam $motherBase.bam $motherBase.fixmate.bam &&
";
print SUB "samtools sort -@ $threads_0 -O bam -o $motherBase.sorted.bam $motherBase.fixmate.bam &&
";
###Remove bam files except *.sorted.bam
print SUB "ls | grep -v sorted | xargs rm
";
=cut

###Child
system "bwa mem -t $threads $chrM_ref $childFull1 $childFull2 -R \'\@RG\\tID:$childBase\\tSM:$childBase\\tPL:ILLUMINA\' | samtools view -@ $threads_0 -O bam -o $childBase.bam"
    and die "\e[31m$childBase Mapping Failed: $!\e[m";
system "samtools fixmate -@ $threads_0 -O bam $childBase.bam $childBase.fixmate.bam"
    and die "\e[31m$childBase Fixmate Failed: $!\e[m";
system "samtools sort -@ $threads_0 -O bam -o $childBase.sorted.bam $childBase.fixmate.bam"
    and die "\e[31m$childBase Sort Failed: $!\e[m";
system "java -jar $picard_path MarkDuplicates I=$childBase.sorted.bam O=$childBase.dedup.bam M=$childBase.metrics.txt"
    and die "\e[31m$childBase mark duplicates Failed: $!\e[m";

###Mother
system "bwa mem -t $threads $chrM_ref $motherFull1 $motherFull2 -R \'\@RG\\tID:$motherBase\\tSM:$motherBase\\tPL:ILLUMINA\' | samtools view -@ $threads_0 -O bam -o $motherBase.bam"
    and die "\e[31m$motherBase Mapping Failed: $!\e[m";
system "samtools fixmate -@ $threads_0 -O bam $motherBase.bam $motherBase.fixmate.bam"
    and die "\e[31m$motherBase Fixmate Failed: $!\e[m";
system "samtools sort -@ $threads_0 -O bam -o $motherBase.sorted.bam $motherBase.fixmate.bam"
    and die "\e[31m$motherBase Sort Failed: $!\e[m";
system "java -jar $picard_path MarkDuplicates I=$motherBase.sorted.bam O=$motherBase.dedup.bam M=$motherBase.metrics.txt"
    and die "\e[31m$motherBase mark duplicates Failed: $!\e[m";

##Delete all bam files except *.dedup.bam
unlink grep {!/dedup/} glob "*.bam";

PILEUP:
##Pileup to generate raw result
chdir($workdir . "1_Mapping") or die "\e[31mCannot enter ${workdir}1_Mapping\e[m";
open my $child_parsed, '>', "${workdir}/2_Pileup/$childBase.parsed";
open my $mother_parsed, '>', "${workdir}/2_Pileup/$motherBase.parsed";
print $child_parsed parse(samtools_mpileup("$childBase.dedup.bam"));
print $mother_parsed parse(samtools_mpileup("$motherBase.dedup.bam"));
close $child_parsed;
close $mother_parsed;

ANNO:
##Annotation
chdir "$workdir/2_Pileup";
open my $childAnno, '>', "$workdir/3_Annotation/$childBase.Annotated.child.tab";
open my $motherAnno, '>', "$workdir/3_Annotation/$motherBase.Annotated.mother.tab";
my ($childAnnoOut, $childPositionDetail) = annotate(homo_het(indel_norm("$childBase.parsed")));
print $childAnno $childAnnoOut;
my ($motherAnnoOut, $motherPositionDetail) = annotate(homo_het(indel_norm("$motherBase.parsed")));
print $motherAnno $motherAnnoOut;
close $childAnno;
close $motherAnno;

##Compare child and mother loci
open my $final_fl, '>', "$workdir/3_Annotation/$childBase.Annotated.merged.tab";
print $final_fl compare("$workdir/3_Annotation/$childBase.Annotated.child.tab",
						"$workdir/3_Annotation/$motherBase.Annotated.mother.tab",
						$childPositionDetail,
						$motherPositionDetail,
					   );


##subroutine for CONFIG readin
sub ReadCONFIG {
    my $line = shift @_;
    $line =~ m/(\S+)\s*=\s*(\S+)/;
    return $2;
}

##Easy to interpret \e[m format to surround text
sub ColorText {
    my ($colorIndex, $text) = @_;
    return("\e[${colorIndex}m$text\e[m");
}
