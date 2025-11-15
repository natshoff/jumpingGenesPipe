#!/usr/bin/env perl -w
use strict;

my @dir;
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	push (@dir, $_);
}
my %values;
my %types;
my %count_for_avg;
my %variance;
for my $dir_s (@dir){
	chdir($dir_s);
my @files = <4_Transposome/*out/*report_annotations_summary.tsv>;
my $total_counts=0;
CHECK: for my $file (@files){
		$file =~ /(.*?)_transposome_(\d+)_out/;
		print "$file\n";
		my $name = $1;
		my $counter=$2;
		my $name_short =$name;
		if($name_short){
			print "$name_short\n";
		}
		else{
			die;
		}
		my %total_values;
		$name=$name."_".$counter;
		open my $tfile, "<", $file;
		while(<$tfile>){
			if(/Species/){
				next;
			}
			if($total_counts == 100){
				last CHECK;
			}
			chomp;
			if(/ReadNum/){
				next;
			}
			my @tarray=split/\s+/;
			
			$values{$name_short}{$tarray[2]}+=$tarray[4];
			push(@{$total_values{$name_short}{$tarray[2]}}, $tarray[4]);
			$types{$tarray[2]}=1;
		}
		$total_counts++;
                        $count_for_avg{$name_short}++;
			for my $typeid (keys %{$total_values{$name_short}}){
				my $temp_total=0;
				for my $valid (@{$total_values{$name_short}{$typeid}}){
					$temp_total+=$valid;
				}
				push (@{$variance{$name_short}{$typeid}}, $temp_total);
			}


}
chdir("../");
}
open my $out, ">", $ARGV[1] . "_TE_annotation_results_average.txt";
open my $out2, ">", $ARGV[1] . "_TE_annotation_results_stdev.txt";
print $out "Species\t";
print $out2 "Species\t";
for my $id(sort keys %types){
	print $out2 "\t$id";
	print $out "\t$id";
}
print $out "\n";
print $out2 "\n";
for my $id (sort keys %values){
		print "$id\n";
		print $out "$id\t";
		print $out2 "$id\t";
		for my $fam (sort keys %types){ 
                if(!exists $values{$id}{$fam}){
			print $out "0\t";
			print $out2 "0\t";
		}
		else{
			my $avg_temp = ($values{$id}{$fam}/$count_for_avg{$id});
			print $out "$avg_temp\t";
			my $tvar;
			for my $varval (@{$variance{$id}{$fam}}){
				$tvar+=($varval-$avg_temp)**2;
			}
			my $var_true = $tvar/$count_for_avg{$id};
			$var_true = sqrt($var_true);
			print $out2 "$var_true\t";
		}
		}
		print $out "\n";
		print $out2 "\n";
}

