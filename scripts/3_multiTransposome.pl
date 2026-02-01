#!/usr/env perl

# NP Hofford Feb 10 2024
# This is the transposome submission script which determines the number of reads sampled and number of Transposome runs ro carry out


use warnings;
use strict;

my $total_pairs = `grep ">" -c $ARGV[0]`;
my $id = $ARGV[1];
# CHANGE numbers here for number of transposome runs
# ex) 100 runs --> (my $i=0; $i<=99; $i++)
for (my $i=0; $i<=99; $i++){
        my $start_rand;
        until($start_rand && $start_rand%4==0){
		# CHANGE number here for the total number of reads sampled
		# ex) $total_pairs-200000
                $start_rand = int(rand($total_pairs-200000));
        }
	# CHANGE head -n XXXXXX number to match total pairs
        `tail -n $start_rand $ARGV[0] | head -n 200000 > $id\_downsample_$start_rand.fasta`;
        `cp ../../../../scripts/GRStransposome_config.yml $id\_$start_rand\_transposome.yml`;
        my $sample_file = $id . "_downsample_$start_rand.fasta";
        my $filePath=$ARGV[0];
        $filePath =~ s/\//\\\//g;
        `perl -pi -e "s/FILE/$sample_file/" $id\_$start_rand\_transposome.yml`;
        my $pwd = `pwd`;
        chomp ($pwd);
        $pwd = $id . "_transposome_$start_rand\_out";
        $pwd =~ s/\//\\\//g;
        `perl -pi -e "s/OUTPUT/$pwd/" $id\_$start_rand\_transposome.yml`;
        `perl -pi -e "s/SAMPLE/$id/" $id\_$start_rand\_transposome.yml`;
        my $config = `pwd`;
        chomp ($config);
        $config .= "/" . $id . "_" . $start_rand . "_transposome.yml";
        $config =~ s/\//\\\//g;
        `sbatch ../../../../scripts/4_transposome.srun $id\_$start_rand\_transposome.yml`;
        }
