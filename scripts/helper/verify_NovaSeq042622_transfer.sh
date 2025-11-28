#!/bin/bash
#SBATCH -J verify_NovaSeq042622_transfer
#SBATCH -n 1
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH --mem=100G
#SBATCH -o logs/slurm-%j.out
#SBATCH -e logs/slurm-%j.err


# Script to verify NovaSeq.042622 samples exist in both vault and scratch locations
# NP Hofford - Nov 2024

SAMPLE_MAP="samples/samples_locations_NovaSeq042622.txt"
VAULT_DIR="/grps2/mrmckain/Sequence_Vault/NovaSeq.042622"
SCRATCH_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"

echo "=== NovaSeq.042622 Sample Transfer Verification ==="
echo ""
echo "Sample map: $SAMPLE_MAP"
echo "Vault location: $VAULT_DIR"
echo "Scratch location: $SCRATCH_DIR"
echo ""

# Check if sample map exists
if [ ! -f "$SAMPLE_MAP" ]; then
    echo "ERROR: Sample map file not found: $SAMPLE_MAP"
    exit 1
fi

# Counters
TOTAL=0
VAULT_FOUND=0
VAULT_MISSING=0
SCRATCH_FOUND=0
SCRATCH_MISSING=0
BOTH_FOUND=0

echo "Checking samples..."
echo ""

while IFS=$'\t' read -r DIR SAMPLE; do
    # Skip empty lines
    [[ -z "$SAMPLE" ]] && continue
    
    TOTAL=$((TOTAL + 1))
    
    # Check vault location (original .bz2 files)
    VAULT_R1=$(find "$VAULT_DIR" -maxdepth 1 -name "*${SAMPLE}*R1*" ! -name "*.html" 2>/dev/null | head -n 1)
    VAULT_R2=$(find "$VAULT_DIR" -maxdepth 1 -name "*${SAMPLE}*R2*" ! -name "*.html" 2>/dev/null | head -n 1)
    
    # Check scratch location (converted .gz files)
    SCRATCH_R1=$(find "$SCRATCH_DIR" -maxdepth 1 -name "*${SAMPLE}*R1*" ! -name "*.html" 2>/dev/null | head -n 1)
    SCRATCH_R2=$(find "$SCRATCH_DIR" -maxdepth 1 -name "*${SAMPLE}*R2*" ! -name "*.html" 2>/dev/null | head -n 1)
    
    VAULT_OK=false
    SCRATCH_OK=false
    
    # Check vault
    if [ -n "$VAULT_R1" ] && [ -n "$VAULT_R2" ]; then
        VAULT_FOUND=$((VAULT_FOUND + 1))
        VAULT_OK=true
    else
        VAULT_MISSING=$((VAULT_MISSING + 1))
        echo "✗ VAULT MISSING: $SAMPLE"
    fi
    
    # Check scratch
    if [ -n "$SCRATCH_R1" ] && [ -n "$SCRATCH_R2" ]; then
        SCRATCH_FOUND=$((SCRATCH_FOUND + 1))
        SCRATCH_OK=true
    else
        SCRATCH_MISSING=$((SCRATCH_MISSING + 1))
        echo "✗ SCRATCH MISSING: $SAMPLE"
    fi
    
    # Check if both found
    if [ "$VAULT_OK" = true ] && [ "$SCRATCH_OK" = true ]; then
        BOTH_FOUND=$((BOTH_FOUND + 1))
        echo "✓ $SAMPLE"
    fi
    
done < "$SAMPLE_MAP"

echo ""
echo "=== Summary ==="
echo "Total samples: $TOTAL"
echo ""
echo "Vault ($VAULT_DIR):"
echo "  Found: $VAULT_FOUND"
echo "  Missing: $VAULT_MISSING"
echo ""
echo "Scratch ($SCRATCH_DIR):"
echo "  Found: $SCRATCH_FOUND"
echo "  Missing: $SCRATCH_MISSING"
echo ""
echo "Both locations OK: $BOTH_FOUND / $TOTAL"
echo ""

if [ $BOTH_FOUND -eq $TOTAL ]; then
    echo "✓ SUCCESS: All samples verified in both locations!"
    exit 0
else
    echo "✗ WARNING: Some samples are missing. See above for details."
    exit 1
fi

