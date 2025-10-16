#!/bin/bash

# ETT_score_2.sh - PRS Scoring Script
# This processes all PRS score files for a given sample
# Usage: ./ETT_score_2.sh <input_prefix> <output_dir> [prs_score_dir]

# Parse arguments
INPUT_PREFIX=$1
OUTPUT_DIR=$2
PRS_SCORE_DIR=${3:-"/Users/alperbulbul/Desktop/PRS"}  # Default PRS directory

# Check arguments
if [ -z "$INPUT_PREFIX" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_prefix> <output_dir> [prs_score_dir]"
    echo ""
    echo "Arguments:"
    echo "  input_prefix  : Path prefix for PLINK files (without .bed/.bim/.fam)"
    echo "  output_dir    : Directory where score results will be saved"
    echo "  prs_score_dir : Directory containing .par files (optional)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/HB0402_id /path/to/output /path/to/PRS"
    exit 1
fi

# PLINK2 binary path (update if needed)
PLINK2=${PLINK2:-"plink2"}  # Use environment variable or default to "plink2"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if input files exist
if [ ! -f "${INPUT_PREFIX}.bed" ] || [ ! -f "${INPUT_PREFIX}.bim" ] || [ ! -f "${INPUT_PREFIX}.fam" ]; then
    echo -e "${RED}Error: Input PLINK files not found at ${INPUT_PREFIX}${NC}"
    echo "Expected files:"
    echo "  ${INPUT_PREFIX}.bed"
    echo "  ${INPUT_PREFIX}.bim"
    echo "  ${INPUT_PREFIX}.fam"
    exit 1
fi

# Check if PRS score directory exists
if [ ! -d "$PRS_SCORE_DIR" ]; then
    echo -e "${RED}Error: PRS score directory not found: $PRS_SCORE_DIR${NC}"
    exit 1
fi

# Count .par files
PAR_COUNT=$(ls -1 "$PRS_SCORE_DIR"/*.par 2>/dev/null | wc -l)
if [ "$PAR_COUNT" -eq 0 ]; then
    echo -e "${RED}Error: No .par files found in $PRS_SCORE_DIR${NC}"
    exit 1
fi

echo "=========================================="
echo "PRS Pipeline Score Calculation"
echo "=========================================="
echo "Input prefix: $INPUT_PREFIX"
echo "Output directory: $OUTPUT_DIR"
echo "PRS score directory: $PRS_SCORE_DIR"
echo "Found $PAR_COUNT score files to process"
echo "=========================================="
echo ""

# Counter for successful/failed runs
SUCCESS_COUNT=0
FAIL_COUNT=0
FAILED_SCORES=""

# Process each score file
for PAR_FILE in "$PRS_SCORE_DIR"/*.par; do
    # Get the base name without extension
    SCORE_NAME=$(basename "$PAR_FILE" .par)
    
    echo -e "${YELLOW}Processing: $SCORE_NAME${NC}"
    
    # Run PLINK2 scoring
    $PLINK2 \
        --bfile "$INPUT_PREFIX" \
        --out "$OUTPUT_DIR/${SCORE_NAME}" \
        --score "$PAR_FILE" 1 2 4 \
        header \
        cols=+scoresums \
        no-mean-imputation \
        2>&1 | tee "$OUTPUT_DIR/${SCORE_NAME}.log"
    
    # Check if successful
    if [ $? -eq 0 ] && [ -f "$OUTPUT_DIR/${SCORE_NAME}.sscore" ]; then
        echo -e "${GREEN}✓ $SCORE_NAME completed successfully${NC}"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        
        # Show brief summary of the score
        if [ -f "$OUTPUT_DIR/${SCORE_NAME}.sscore" ]; then
            LINES=$(wc -l < "$OUTPUT_DIR/${SCORE_NAME}.sscore")
            echo "  Output: ${SCORE_NAME}.sscore ($LINES lines)"
        fi
    else
        echo -e "${RED}✗ $SCORE_NAME failed${NC}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        FAILED_SCORES="$FAILED_SCORES $SCORE_NAME"
    fi
    
    echo "----------------------------------------"
done

# Summary report
echo ""
echo "=========================================="
echo "SUMMARY REPORT"
echo "=========================================="
echo -e "${GREEN}Successful: $SUCCESS_COUNT scores${NC}"
echo -e "${RED}Failed: $FAIL_COUNT scores${NC}"

if [ $FAIL_COUNT -gt 0 ]; then
    echo ""
    echo "Failed scores:$FAILED_SCORES"
fi

echo ""
echo "Results saved in: $OUTPUT_DIR"
echo ""

# List all generated score files
echo "Generated score files:"
ls -lh "$OUTPUT_DIR"/*.sscore 2>/dev/null | awk '{print "  - " $9 " (" $5 ")"}'

echo ""
echo "=========================================="

# Exit with appropriate code
if [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${GREEN}All scores processed successfully!${NC}"
    exit 0
else
    echo -e "${YELLOW}Processing completed with some failures${NC}"
    exit 1
fi