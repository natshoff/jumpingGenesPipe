#!/bin/bash

# Quick script to list what samples are actually in the scratch directory

SCRATCH_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"

echo "=== Samples in Scratch Directory ==="
echo ""

# Find all R1 files and extract sample IDs
ls -1 "$SCRATCH_DIR"/*R1*.gz 2>/dev/null | while read file; do
    basename "$file" | sed 's/.*JG_CO_/JG_CO_/' | sed 's/_I[0-9].*//g' | sed 's/_L[0-9].*//g'
done | sort -u

echo ""
echo "Total unique samples:"
ls -1 "$SCRATCH_DIR"/*R1*.gz 2>/dev/null | wc -l

