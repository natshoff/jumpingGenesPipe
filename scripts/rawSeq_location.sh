#!/bin/bash
#SBATCH -J sample_map
#SBATCH -n 1
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH --mem=4G
#SBATCH -t 5:00:00
#SBATCH -o logs/rawSeq_loc-%j.out
#SBATCH -e logs/rawSeq_loc-%j.err

# Script to generate a sample location map
# Creates a two-column file: directory_path | sample_id
# 
# Usage: sbatch rawSeq_loc.sh samples/samples.txt samples/sample_locations.txt
# Or run directly: ./rawSeq_loc.sh samples/samples.txt samples/sample_locations.txt

if [ $# -lt 2 ]; then
    echo "Usage: $0 SAMPLE_LIST OUTPUT_FILE"
    echo ""
    echo "Example: $0 samples/samples.txt samples/sample_locations.txt"
    exit 1
fi

SAMPLE_LIST=$1
OUTPUT_FILE=$2

# Create logs directory if needed
mkdir -p logs

# Define sequence directories to search
# Note: NovaSeq.042622 has .bz2 files that are converted to .gz in scratch
SEARCH_DIRS=(
    "/grps2/mrmckain/Sequence_Vault/NovaSeq.022122"
    "/grps2/mrmckain/Sequence_Vault/NovaSeq.042622"
    "/grps2/mrmckain/Sequence_Vault/NovaSeq.110321"
)

echo "Generating sample location map..."
echo "Input: $SAMPLE_LIST"
echo "Output: $OUTPUT_FILE"
echo ""

# Clear output file if it exists
> "$OUTPUT_FILE"

# Track statistics
FOUND=0
NOT_FOUND=0

# Read each sample from the input file
while IFS= read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue
    
    # Clean up whitespace
    SAMPLE=$(echo "$SAMPLE" | xargs)
    
    FOUND_IN=""
    
    # Search each directory
    for DIR in "${SEARCH_DIRS[@]}"; do
        # Check if directory exists
        if [ ! -d "$DIR" ]; then
            continue
        fi
        
        # Look for R1 file matching this sample (exclude .html files)
        R1_FILE=$(find "$DIR" -maxdepth 1 -name "*${SAMPLE}*R1*" ! -name "*.html" 2>/dev/null | head -n 1)
        
        if [ -n "$R1_FILE" ]; then
            FOUND_IN="$DIR"
            break
        fi
    done
    
    if [ -n "$FOUND_IN" ]; then
        # Remap NovaSeq.042622 vault location to scratch location (where converted files are)
        OUTPUT_DIR="$FOUND_IN"
        if [ "$FOUND_IN" == "/grps2/mrmckain/Sequence_Vault/NovaSeq.042622" ]; then
            OUTPUT_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"
            echo "✓ $SAMPLE -> $FOUND_IN (mapped to scratch)"
        else
            echo "✓ $SAMPLE -> $FOUND_IN"
        fi
        
        echo "$OUTPUT_DIR	$SAMPLE" >> "$OUTPUT_FILE"
        FOUND=$((FOUND + 1))
    else
        echo "✗ $SAMPLE -> NOT FOUND"
        NOT_FOUND=$((NOT_FOUND + 1))
    fi
    
done < "$SAMPLE_LIST"

echo ""
echo "=== Summary ==="
echo "Found: $FOUND samples"
echo "Not found: $NOT_FOUND samples"
echo ""
echo "Sample location map written to: $OUTPUT_FILE"

