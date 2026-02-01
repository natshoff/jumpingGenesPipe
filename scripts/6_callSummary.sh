#!/bin/bash
#SBATCH -n 4
#SBATCH -p main
#SBATCH -N 1
#SBATCH --mem=25G
#SBATCH --job-name=TE_summary
#SBATCH --output=TE_summary_%j.out
#SBATCH --error=TE_summary_%j.err

# 6_callSummary.sh
# Calls 7_pullAnnotation.pl to generate TE annotation summaries
#
# Usage: sbatch 6_callSummary.sh <data_path_file> <level> <output_prefix>
#
# Arguments:
#   data_path_file  - Text file containing paths to sample directories (from 5_processedDataPath.pl)
#   level           - Classification level: "order", "superfamily", or "family"
#   output_prefix   - Prefix for output files (will generate _average.txt and _stdev.txt)
#
# Example:
#   sbatch 6_callSummary.sh samples.txt superfamily my_analysis
#
# Output files:
#   <output_prefix>_<level>_TE_annotation_results_average.txt
#   <output_prefix>_<level>_TE_annotation_results_stdev.txt

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: sbatch 6_callSummary.sh <data_path_file> <level> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  data_path_file  - Text file with sample directory paths (from 5_processedDataPath.pl)"
    echo "  level           - Classification level: order, superfamily, or family"
    echo "  output_prefix   - Prefix for output files"
    echo ""
    echo "Example:"
    echo "  sbatch 6_callSummary.sh samples.txt superfamily sorghum_analysis"
    exit 1
fi

DATA_PATH_FILE=$1
LEVEL=$2
OUTPUT_PREFIX=$3

# Validate level argument
if [[ "$LEVEL" != "order" && "$LEVEL" != "superfamily" && "$LEVEL" != "family" ]]; then
    echo "Error: level must be 'order', 'superfamily', or 'family'"
    echo "You provided: $LEVEL"
    exit 1
fi

echo "=========================================="
echo "TE Annotation Summary Pipeline"
echo "=========================================="
echo "Data path file: $DATA_PATH_FILE"
echo "Classification level: $LEVEL"
echo "Output prefix: $OUTPUT_PREFIX"
echo "=========================================="
echo ""

# Run the Perl script
perl "${SCRIPT_DIR}/7_pullAnnotation.pl" "$DATA_PATH_FILE" "$LEVEL" "$OUTPUT_PREFIX"

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_${LEVEL}_TE_annotation_results_average.txt"
echo "  ${OUTPUT_PREFIX}_${LEVEL}_TE_annotation_results_stdev.txt"
echo "=========================================="
