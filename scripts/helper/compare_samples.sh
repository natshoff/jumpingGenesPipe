#!/bin/bash

# Compare samples in location map vs actual scratch directory

SAMPLE_MAP="samples/samples_locations_NovaSeq042622.txt"
SCRATCH_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"

echo "=== Sample Comparison ==="
echo ""

# Extract sample IDs from location map
echo "Samples in location map:"
awk '{print $2}' "$SAMPLE_MAP" | grep -v "^$" | sort > /tmp/map_samples.txt
MAP_COUNT=$(cat /tmp/map_samples.txt | wc -l)
echo "  Count: $MAP_COUNT"
echo ""

# Extract sample IDs from scratch directory (from R1 files)
echo "Samples in scratch directory:"
ls -1 "$SCRATCH_DIR"/*R1* 2>/dev/null | while read file; do
    basename "$file" | sed 's/.*-\(JG_CO_[^_]*_[^_]*\).*/\1/' | sed 's/.*_HA_//' | sed 's/_I[0-9].*//'
done | sort -u > /tmp/scratch_samples.txt
SCRATCH_COUNT=$(cat /tmp/scratch_samples.txt | wc -l)
echo "  Count: $SCRATCH_COUNT"
echo ""

# Find samples in map but NOT in scratch
echo "=== Samples in map but MISSING from scratch ==="
comm -23 /tmp/map_samples.txt /tmp/scratch_samples.txt
MISSING_COUNT=$(comm -23 /tmp/map_samples.txt /tmp/scratch_samples.txt | wc -l)
echo "  Count: $MISSING_COUNT"
echo ""

# Find samples in scratch but NOT in map
echo "=== Samples in scratch but NOT in map ==="
comm -13 /tmp/map_samples.txt /tmp/scratch_samples.txt
EXTRA_COUNT=$(comm -13 /tmp/map_samples.txt /tmp/scratch_samples.txt | wc -l)
echo "  Count: $EXTRA_COUNT"
echo ""

echo "=== Summary ==="
echo "Location map: $MAP_COUNT samples"
echo "Scratch dir:  $SCRATCH_COUNT samples"
echo "Missing:      $MISSING_COUNT samples"
echo "Extra:        $EXTRA_COUNT samples"

# Cleanup
rm -f /tmp/map_samples.txt /tmp/scratch_samples.txt

