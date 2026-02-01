#!/usr/bin/env perl -w
# 5_processedDataPath.pl
# Generates a list of valid sample directories for the annotation pipeline
#
# Usage: perl 5_processedDataPath.pl <base_directory> <output_prefix> <expected_count>
#
# Arguments:
#   base_directory  - Path to processed data (e.g., /scratch/.../data/processed)
#   output_prefix   - Prefix for output files
#   expected_count  - Expected number of transposome output files (e.g., 100)
#
# Example:
#   perl 5_processedDataPath.pl /scratch/nphofford/jumpingGenesPipe/data/processed samples 100
#
# Output files:
#   <output_prefix>_valid.txt    - Samples with exactly <expected_count> annotation files
#   <output_prefix>_excluded.txt - Samples with different counts (includes count in file)

use strict;
use File::Find;
use File::Basename;

# Check arguments
if (@ARGV < 3) {
    print STDERR "Usage: perl 5_processedDataPath.pl <base_directory> <output_prefix> <expected_count>\n";
    print STDERR "\n";
    print STDERR "Arguments:\n";
    print STDERR "  base_directory  - Path to processed data directory\n";
    print STDERR "  output_prefix   - Prefix for output files\n";
    print STDERR "  expected_count  - Expected number of transposome output files (e.g., 100)\n";
    print STDERR "\n";
    print STDERR "Example:\n";
    print STDERR "  perl 5_processedDataPath.pl /scratch/nphofford/jumpingGenesPipe/data/processed samples 100\n";
    print STDERR "\n";
    print STDERR "Output:\n";
    print STDERR "  samples_valid.txt    - Paths with exactly 100 annotation files\n";
    print STDERR "  samples_excluded.txt - Paths with different counts\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $output_prefix = $ARGV[1];
my $expected_count = $ARGV[2];

# Validate expected_count is a positive integer
unless ($expected_count =~ /^\d+$/ && $expected_count > 0) {
    die "Error: expected_count must be a positive integer. You provided: $expected_count\n";
}

my $valid_file = "${output_prefix}_valid.txt";
my $excluded_file = "${output_prefix}_excluded.txt";

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
my $excluded_count = 0;
my $no_transposome_count = 0;
my @valid_samples;
my @excluded_samples;

print "Scanning for sample directories in: $base_dir\n";
print "Expected annotation file count: $expected_count\n";
print "=" x 70 . "\n";

for my $sample (@sample_dirs) {
    $total_count++;
    my $sample_path = "$base_dir/$sample";
    my $transposome_dir = "$sample_path/4_Transposome";
    
    # Check if 4_Transposome directory exists
    unless (-d $transposome_dir) {
        print "[SKIP] $sample - no 4_Transposome directory\n";
        push @excluded_samples, { path => $sample_path, count => 0, reason => "no_transposome_dir" };
        $no_transposome_count++;
        next;
    }
    
    # Find annotation summary files
    my @tsv_files = glob("$transposome_dir/*out/*report_annotations_summary.tsv");
    my $tsv_count = scalar @tsv_files;
    
    if ($tsv_count == $expected_count) {
        # Exactly the expected count - valid
        push @valid_samples, $sample_path;
        $valid_count++;
        print "[OK]   $sample ($tsv_count files)\n";
    } elsif ($tsv_count == 0) {
        # No files at all
        print "[EXCL] $sample - no annotation files\n";
        push @excluded_samples, { path => $sample_path, count => $tsv_count, reason => "no_files" };
        $excluded_count++;
    } else {
        # Has files but not the expected count
        my $diff = $tsv_count - $expected_count;
        my $diff_str = $diff > 0 ? "+$diff" : "$diff";
        print "[EXCL] $sample ($tsv_count files, $diff_str from expected)\n";
        push @excluded_samples, { path => $sample_path, count => $tsv_count, reason => "wrong_count" };
        $excluded_count++;
    }
}

# Write valid samples file
open my $out, ">", $valid_file or die "Cannot open output file $valid_file: $!";
for my $path (@valid_samples) {
    print $out "$path\n";
}
close $out;

# Write excluded samples file (with counts for debugging)
open my $out2, ">", $excluded_file or die "Cannot open output file $excluded_file: $!";
print $out2 "# Samples excluded from analysis\n";
print $out2 "# Expected count: $expected_count\n";
print $out2 "# Format: count<TAB>reason<TAB>path\n";
print $out2 "#\n";
for my $sample (@excluded_samples) {
    print $out2 "$sample->{count}\t$sample->{reason}\t$sample->{path}\n";
}
close $out2;

# Summary
print "=" x 70 . "\n";
print "Summary:\n";
print "  Total directories scanned:    $total_count\n";
print "  Valid (exactly $expected_count files):   $valid_count\n";
print "  Excluded (wrong count):       $excluded_count\n";
print "  No 4_Transposome directory:   $no_transposome_count\n";
print "\n";
print "Output files:\n";
print "  Valid samples:    $valid_file\n";
print "  Excluded samples: $excluded_file\n";
print "\n";
print "Next step:\n";
print "  sbatch 6_callSummary.sh $valid_file <order|superfamily|family> <output_prefix>\n";
