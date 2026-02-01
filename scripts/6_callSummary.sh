#!/bin/bash
#SBATCH -n 4
#SBATCH --qos mrmckain
#SBATCH -p ultrahigh
#SBATCH -N 1
#SBATCH --mem=25G
#SBATCH --job-name=TE_summary
#SBATCH --output=TE_summary_%j.out
#SBATCH --error=TE_summary_%j.err

# 6_callSummary.sh
# Calls 7_pullAnnotation.pl to generate TE annotation summaries
#
# Usage: sbatch 6_callSummary.sh <data_path_file> <level> <output_name>
#
# Arguments:
#   data_path_file  - Text file containing paths to sample directories (from 5_processedDataPath.pl)
#   level           - Classification level: "order", "superfamily", or "family"
#   output_name     - Name for output files (e.g., feb1superfamily)
#
# Example:
#   sbatch 6_callSummary.sh /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt superfamily feb1superfamily
#
# Output files (in /scratch/nphofford/jumpingGenesPipe/data/output/):
#   <output_name>_<level>_TE_annotation_results_average.txt
#   <output_name>_<level>_TE_annotation_results_stdev.txt

# Hardcoded paths
SCRIPT_DIR="/scratch/nphofford/jumpingGenesPipe/scripts"
OUTPUT_DIR="/scratch/nphofford/jumpingGenesPipe/data/output"

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: sbatch 6_callSummary.sh <data_path_file> <level> <output_name>"
    echo ""
    echo "Arguments:"
    echo "  data_path_file  - Text file with sample directory paths (from 5_processedDataPath.pl)"
    echo "  level           - Classification level: order, superfamily, or family"
    echo "  output_name     - Name for output files (e.g., feb1superfamily)"
    echo ""
    echo "Example:"
    echo "  sbatch 6_callSummary.sh /scratch/nphofford/jumpingGenesPipe/samples/outputRuns/feb1samples_valid.txt superfamily feb1superfamily"
    echo ""
    echo "Output location: $OUTPUT_DIR/"
    exit 1
fi

DATA_PATH_FILE=$1
LEVEL=$2
OUTPUT_NAME=$3

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
echo "Output name: $OUTPUT_NAME"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="
echo ""

# Run the Perl script
perl "${SCRIPT_DIR}/7_pullAnnotation.pl" "$DATA_PATH_FILE" "$LEVEL" "$OUTPUT_NAME"

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "Output files:"
echo "  ${OUTPUT_DIR}/${OUTPUT_NAME}_${LEVEL}_TE_annotation_results_average.txt"
echo "  ${OUTPUT_DIR}/${OUTPUT_NAME}_${LEVEL}_TE_annotation_results_stdev.txt"
echo "=========================================="
