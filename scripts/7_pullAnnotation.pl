#!/usr/bin/env perl -w
# 7_pullAnnotation.pl
# Pulls and summarizes TE annotations from transposome output
# Calculates average and standard deviation of GenomeFraction across runs
#
# Usage: perl 7_pullAnnotation.pl <data_path_file> <level> <output_prefix>
#
# Arguments:
#   data_path_file  - Text file containing paths to sample directories
#   level           - Classification level: "order", "superfamily", or "family"
#   output_prefix   - Prefix for output files
#
# Output:
#   <output_prefix>_<level>_TE_annotation_results_average.txt
#   <output_prefix>_<level>_TE_annotation_results_stdev.txt
#
# TSV Column Reference:
#   0: ReadNum
#   1: Order (e.g., ltr_retrotransposon, dna_transposon, satellite)
#   2: Superfamily (e.g., Gypsy, Copia)
#   3: Family (e.g., ATHILA-1_SBi-I, Copia-141_SB-I)
#   4: GenomeFraction

use strict;

# Check arguments
if (@ARGV < 3) {
    print STDERR "Usage: perl 7_pullAnnotation.pl <data_path_file> <level> <output_prefix>\n";
    print STDERR "\n";
    print STDERR "Arguments:\n";
    print STDERR "  data_path_file  - Text file with sample directory paths\n";
    print STDERR "  level           - Classification level: order, superfamily, or family\n";
    print STDERR "  output_prefix   - Prefix for output files\n";
    print STDERR "\n";
    print STDERR "Example:\n";
    print STDERR "  perl 7_pullAnnotation.pl samples.txt superfamily sorghum_analysis\n";
    exit 1;
}

my $data_path_file = $ARGV[0];
my $level = lc($ARGV[1]);  # Normalize to lowercase
my $output_prefix = $ARGV[2];

# Determine which column to use based on level
my $level_column;
if ($level eq "order") {
    $level_column = 1;
} elsif ($level eq "superfamily") {
    $level_column = 2;
} elsif ($level eq "family") {
    $level_column = 3;
} else {
    die "Error: level must be 'order', 'superfamily', or 'family'. You provided: $level\n";
}

print "Classification level: $level (column $level_column)\n";
print "GenomeFraction column: 4\n\n";

# Read sample directories from file
my @dir;
open my $file, "<", $data_path_file or die "Cannot open $data_path_file: $!";
while (<$file>) {
    chomp;
    next if /^\s*$/;  # Skip empty lines
    push(@dir, $_);
}
close $file;

print "Loaded " . scalar(@dir) . " sample directories\n\n";

# Data structures for aggregation
my %values;         # Sum of GenomeFraction per sample per type
my %types;          # All unique type names encountered
my %count_for_avg;  # Number of runs per sample
my %variance;       # Values for standard deviation calculation

# Process each sample directory
for my $dir_s (@dir) {
    chdir($dir_s) or do {
        warn "Warning: Cannot change to directory $dir_s: $!\n";
        next;
    };
    
    my @files = <4_Transposome/*out/*report_annotations_summary.tsv>;
    
    if (@files == 0) {
        warn "Warning: No annotation files found in $dir_s/4_Transposome/\n";
        chdir("../");
        next;
    }
    
    print "Processing: $dir_s (" . scalar(@files) . " runs)\n";
    
    for my $file (@files) {
        # Extract sample name from file path
        # Pattern: SAMPLENAME_transposome_JOBID_out
        $file =~ /(.*?)_transposome_(\d+)_out/;
        my $name = $1;
        my $counter = $2;
        
        my $name_short = $name;
        unless ($name_short) {
            warn "Warning: Could not parse sample name from $file\n";
            next;
        }
        
        my %total_values;  # Temporary storage for this run
        my $name_full = "${name}_${counter}";
        
        open my $tfile, "<", $file or do {
            warn "Warning: Cannot open $file: $!\n";
            next;
        };
        
        while (<$tfile>) {
            # Skip header lines
            next if /Species/;
            next if /ReadNum/;
            
            chomp;
            my @tarray = split /\s+/;
            
            # Validate we have enough columns
            next unless defined $tarray[4];
            
            my $type_name = $tarray[$level_column];
            my $genome_fraction = $tarray[4];
            
            # Accumulate values
            $values{$name_short}{$type_name} += $genome_fraction;
            push(@{$total_values{$name_short}{$type_name}}, $genome_fraction);
            $types{$type_name} = 1;
        }
        close $tfile;
        
        # Count this run for the sample
        $count_for_avg{$name_short}++;
        
        # Store totals per run for variance calculation
        for my $typeid (keys %{$total_values{$name_short}}) {
            my $temp_total = 0;
            for my $valid (@{$total_values{$name_short}{$typeid}}) {
                $temp_total += $valid;
            }
            push(@{$variance{$name_short}{$typeid}}, $temp_total);
        }
    }
    
    chdir("../");
}

# Generate output file names
my $avg_file = "${output_prefix}_${level}_TE_annotation_results_average.txt";
my $stdev_file = "${output_prefix}_${level}_TE_annotation_results_stdev.txt";

print "\nWriting output files...\n";

# Open output files
open my $out, ">", $avg_file or die "Cannot open $avg_file: $!";
open my $out2, ">", $stdev_file or die "Cannot open $stdev_file: $!";

# Write headers
print $out "Species";
print $out2 "Species";
for my $id (sort keys %types) {
    print $out "\t$id";
    print $out2 "\t$id";
}
print $out "\n";
print $out2 "\n";

# Write data rows
for my $id (sort keys %values) {
    print "  Sample: $id ($count_for_avg{$id} runs)\n";
    print $out "$id";
    print $out2 "$id";
    
    for my $fam (sort keys %types) {
        if (!exists $values{$id}{$fam}) {
            print $out "\t0";
            print $out2 "\t0";
        } else {
            # Calculate average
            my $avg_temp = $values{$id}{$fam} / $count_for_avg{$id};
            print $out "\t$avg_temp";
            
            # Calculate standard deviation
            my $tvar = 0;
            for my $varval (@{$variance{$id}{$fam}}) {
                $tvar += ($varval - $avg_temp) ** 2;
            }
            my $var_true = sqrt($tvar / $count_for_avg{$id});
            print $out2 "\t$var_true";
        }
    }
    
    print $out "\n";
    print $out2 "\n";
}

close $out;
close $out2;

print "\n";
print "=" x 50 . "\n";
print "Summary complete!\n";
print "=" x 50 . "\n";
print "Samples processed: " . scalar(keys %values) . "\n";
print "Unique ${level}s found: " . scalar(keys %types) . "\n";
print "\n";
print "Output files:\n";
print "  Average: $avg_file\n";
print "  Std Dev: $stdev_file\n";
