#!/usr/bin/env perl -w
use strict;

$ARGV[0] =~ /(.+).fq/;
my %seq;
my %seq2;
open my $out, ">", $1 . ".fasta";

open my $file, "<", $ARGV[0];
open my $file2, "<", $ARGV[1];
my $count=0;
while(<$file>){
	chomp;
	my $id = $_;
	my $seq=readline($file);
	readline($file);
	readline($file);
	$seq{$count}=">".$id."\n".$seq;
$count++;
	#print $out ">$id\n$seq";
}

$count=0;
while(<$file2>){
        chomp;
        my $id = $_;
        my $seq=readline($file2);
        readline($file2);
        readline($file2);
        $seq2{$count}=">".$id."\n".$seq;
$count++;
        #print $out ">$id\n$seq";
        }

for my $cou (sort {$a<=>$b} keys %seq){
	print $out "$seq{$cou}";
	print $out "$seq2{$cou}";
}
