#!/usr/bin/perl -w

# NP Hofford - Nov 2024
# This script is used to initiate the TE pipeline using a sample location map file
# 
# It is called using three inputs (1) a sample location map file (directory\tsample_id format), (2) the TE_pipeline.sh file, and (3) optional starting step
# example call: perl 1_callTEpipe_mapped.pl samples/sample_locations.txt scripts/2_TEpipe_fastp.sh
# example call: perl 1_callTEpipe_mapped.pl samples/sample_locations.txt scripts/2_TEpipe_fastp.sh 3


use strict;

# Get starting step (default to 1 if not provided)
my $start_step = $ARGV[2] || 1;
print "Starting pipeline from Step $start_step\n\n";

open my $file, "<", $ARGV[0];

while(<$file>){
	chomp;
	my @tarray=split/\t/;
	# Location map format: directory\tsample_id
	my @files = grep { !/\.html$/ } <$tarray[0]/*$tarray[1]*R[1-2]*>;
        print "$tarray[1]\n";
	my $file1;
	my $file2;
	for my $tfile (@files){
		if($tfile =~ /R1/){
			$file1=$tfile;
			print "$file1\n";
		}
		else{
			$file2=$tfile;
			print "$file2\n";
		}
	};
	`sbatch $ARGV[1] $file1 $file2 $tarray[1] $start_step`;
	print "Submitted job for $tarray[1] (starting from Step $start_step)\n\n";
#	chdir("../");
}	
		


