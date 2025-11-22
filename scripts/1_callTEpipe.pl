#!/usr/bin/perl -w

# NP Hofford - Feb 10 2024
# This script is used to initiate the TE pipeline
# 
# It is called using three inputs (1) the path to the raw sequence files (my @files = ), (2) a text file of the samples you want to run, and (3) the TE_pipeline.sh file
# example call: perl samples.txt TE_pipeline.sh 

# Raw sequence files are spread across three sequence runs:
# /grps2/mrmckain/Sequence_Vault/NovaSeq.022122/
# /grps2/mrmckain/Sequence_Vault/NovaSeq.042622/
# /grps2/mrmckain/Sequence_Vault/NovaSeq.110321/


use strict;

open my $file, "<", $ARGV[0];

while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	# CHNAGE YOUR SEQUENCE PATH HERE
	my @files = grep { !/\.html$/ } </scratch/nphofford/Sorghum_halepense_seq/*$tarray[0]*R[1-2].fq*>;
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
	`sbatch $ARGV[1] $file1 $file2 $tarray[0] `;
#	chdir("../");
}	
		
