#!/usr/bin/env perl -w
# 5_processedDataPath.pl
# Generates a list of valid sample directories for the annotation pipeline
#
# Usage: perl 5_processedDataPath.pl <base_directory> <output_file>
#
# Example:
#   perl 5_processedDataPath.pl /scratch/nphofford/jumpingGenesPipe/data/processed samples.txt
#
# Only includes samples that have valid 4_Transposome output with annotation summary files

use strict;
use File::Find;
use File::Basename;

# Check arguments
if (@ARGV < 2) {
    print STDERR "Usage: perl 5_processedDataPath.pl <base_directory> <output_file>\n";
    print STDERR "Example: perl 5_processedDataPath.pl /scratch/nphofford/jumpingGenesPipe/data/processed samples.txt\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $output_file = $ARGV[1];

# Verify base directory exists
unless (-d $base_dir) {
    die "Error: Base directory does not exist: $base_dir\n";
}

# Get all subdirectories in base_dir
opendir(my $dh, $base_dir) or die "Cannot open directory $base_dir: $!";
my @sample_dirs = grep { -d "$base_dir/$_" && $_ !~ /^\./ } readdir($dh);
closedir($dh);

# Sort sample directories
@sample_dirs = sort @sample_dirs;

my $total_count = 0;
my $valid_count = 0;
my @valid_samples;

print "Scanning for valid sample directories in: $base_dir\n";
print "=" x 60 . "\n";

for my $sample (@sample_dirs) {
    $total_count++;
    my $sample_path = "$base_dir/$sample";
    my $transposome_dir = "$sample_path/4_Transposome";
    
    # Check if 4_Transposome directory exists
    unless (-d $transposome_dir) {
        print "[SKIP] $sample - no 4_Transposome directory\n";
        next;
    }
    
    # Find annotation summary files
    my @tsv_files = glob("$transposome_dir/*out/*report_annotations_summary.tsv");
    my $tsv_count = scalar @tsv_files;
    
    if ($tsv_count > 0) {
        push @valid_samples, $sample_path;
        $valid_count++;
        print "[OK]   $sample ($tsv_count annotation files)\n";
    } else {
        print "[SKIP] $sample - no annotation summary files in 4_Transposome\n";
    }
}

# Write output file
open my $out, ">", $output_file or die "Cannot open output file $output_file: $!";
for my $path (@valid_samples) {
    print $out "$path\n";
}
close $out;

# Summary
print "=" x 60 . "\n";
print "Total directories scanned: $total_count\n";
print "Valid samples written: $valid_count\n";
print "Output saved to: $output_file\n";
print "\n";
print "Next step:\n";
print "  sbatch 6_callSummary.sh $output_file <order|superfamily|family> <output_prefix>\n";
