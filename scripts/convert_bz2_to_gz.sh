#!/bin/bash
#SBATCH -J bz2_to_gz
#SBATCH -n 1
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH --mem=20G
#SBATCH -t 10:00:00
#SBATCH -o logs/convert_bz2_to_gz-%j.out
#SBATCH -e logs/convert_bz2_to_gz-%j.err

# Script to convert .bz2 compressed FASTQ files to .gz and copy to scratch directory
# NP Hofford - Nov 2024

# Source directory (only NovaSeq.042622 has .bz2 files)
SOURCE_DIR="/grps2/mrmckain/Sequence_Vault/NovaSeq.042622"

# Destination directory
DEST_DIR="/scratch/nphofford/jumpingGenesPipe/data/raw/NovaSeq.042622"

# Create logs directory if needed
mkdir -p logs

echo "=== Starting .bz2 to .gz conversion ==="
echo "Source: $SOURCE_DIR"
echo "Destination: $DEST_DIR"
echo ""

# Create destination directory
mkdir -p "$DEST_DIR"

# Find all .bz2 files that match JG_CO_### pattern (digits only after JG_CO_)
# Exclude files with H-HiC in the name
echo "Finding .bz2 files with JG_CO_### pattern (excluding H-HiC)..."
BZ2_FILES=$(find "$SOURCE_DIR" -maxdepth 1 \( -name "*JG_CO_[0-9]*.fastq.bz2" -o -name "*JG_CO_[0-9]*.fq.bz2" \) | grep -v "H-HiC")

if [ -z "$BZ2_FILES" ]; then
    echo "No .bz2 files with JG_CO_### pattern found (after filtering)!"
    exit 1
fi

# Count total files
TOTAL=$(echo "$BZ2_FILES" | wc -l)
echo "Found $TOTAL .bz2 files to convert"
echo ""

COUNTER=0
# Convert each .bz2 file
for BZ2_FILE in $BZ2_FILES; do
    COUNTER=$((COUNTER + 1))
    BASENAME=$(basename "$BZ2_FILE" .bz2)
    DEST_FILE="$DEST_DIR/${BASENAME}.gz"
    
    # Skip if already exists
    if [ -f "$DEST_FILE" ]; then
        echo "[$COUNTER/$TOTAL] SKIP: $BASENAME.gz (already exists)"
        continue
    fi
    
    echo "[$COUNTER/$TOTAL] Converting: $(basename "$BZ2_FILE") -> ${BASENAME}.gz"
    
    # Decompress .bz2 and recompress as .gz
    bunzip2 -c "$BZ2_FILE" | gzip -c > "$DEST_FILE"
    
    # Check if conversion was successful
    if [ $? -eq 0 ]; then
        echo "  ✓ Success"
    else
        echo "  ✗ Failed"
        rm -f "$DEST_FILE"  # Remove partial file if failed
    fi
done

echo ""
echo "=== Conversion complete ==="
COUNT=$(find "$DEST_DIR" -name "*.gz" 2>/dev/null | wc -l)
echo "Total .gz files in destination: $COUNT"

