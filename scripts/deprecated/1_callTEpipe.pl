#!/usr/bin/perl -w

# NP Hofford - Feb 10 2024
# This script is used to initiate the TE pipeline
# 
# Arguments:
#   $ARGV[0] - Sample list file
#   $ARGV[1] - Pipeline script (2_TEpipe_fastp.sh)
#   $ARGV[2] - Starting step (optional, default = 1)
#
# Example calls:
#   perl 1_callTEpipe.pl samples.txt scripts/2_TEpipe_fastp.sh
#   perl 1_callTEpipe.pl samples.txt scripts/2_TEpipe_fastp.sh 3

# Raw sequence files are spread across three sequence runs:
# /grps2/mrmckain/Sequence_Vault/NovaSeq.022122/
# /grps2/mrmckain/Sequence_Vault/NovaSeq.042622/
# /grps2/mrmckain/Sequence_Vault/NovaSeq.110321/


use strict;

# Get starting step (default to 1 if not provided)
my $start_step = $ARGV[2] || 1;
print "Starting pipeline from Step $start_step\n\n";

open my $file, "<", $ARGV[0];

while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	# CHNAGE YOUR SEQUENCE PATH HERE
	my @files = grep { !/\.html$/ } </grps2/mrmckain/Sequence_Vault/NovaSeq.022122/*$tarray[0]*R[1-2]*>;
        print "$tarray[0]\n";
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
	`sbatch $ARGV[1] $file1 $file2 $tarray[0] $start_step`;
	print "Submitted job for $tarray[0] (starting from Step $start_step)\n\n";
#	chdir("../");
}	
		
