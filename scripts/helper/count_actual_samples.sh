#!/bin/bash

# Simple script to count actual complete sample pairs in scratch

SCRATCH_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"

echo "=== Actual Sample Count in Scratch ==="
echo ""

# Count total files
TOTAL_FILES=$(ls -1 "$SCRATCH_DIR" 2>/dev/null | wc -l)
echo "Total files: $TOTAL_FILES"

# Count R1 files
R1_COUNT=$(ls -1 "$SCRATCH_DIR"/*R1* 2>/dev/null | wc -l)
echo "R1 files: $R1_COUNT"

# Count R2 files  
R2_COUNT=$(ls -1 "$SCRATCH_DIR"/*R2* 2>/dev/null | wc -l)
echo "R2 files: $R2_COUNT"

echo ""
echo "Expected samples (both R1 and R2): $R1_COUNT"
echo "(Assuming R1 count = R2 count)"

# Extract sample names from R1 files
echo ""
echo "Sample IDs in scratch (from R1 files):"
ls -1 "$SCRATCH_DIR"/*R1* 2>/dev/null | while read f; do
    basename "$f" | sed 's/.*JG_CO_/JG_CO_/' | sed 's/_I[0-9].*//'
done | sort

