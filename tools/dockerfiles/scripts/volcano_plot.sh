#!/bin/bash

# volcano_plot.sh

set -euo pipefail  # Enable strict error handling

#########################
# Argument Parsing      #
#########################

# Initialize variables
DIFF_EXPR_FILE=""
X_AXIS_COLUMN=""
Y_AXIS_COLUMN=""
LABEL_COLUMN=""
OUTPUT_FILENAME="index.html"  # Default value


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input)
            DIFF_EXPR_FILE="$2"
            shift 2
            ;;
        --x)
            X_AXIS_COLUMN="$2"
            shift 2
            ;;
        --y)
            Y_AXIS_COLUMN="$2"
            shift 2
            ;;
        --label)
            LABEL_COLUMN="$2"
            shift 2
            ;;
        --output)
            OUTPUT_FILENAME="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$DIFF_EXPR_FILE" || -z "$X_AXIS_COLUMN" || -z "$Y_AXIS_COLUMN" || -z "$LABEL_COLUMN" ]]; then
    echo "Missing required arguments." >&2
    exit 1
fi

#########################
# Main Script Execution #
#########################

# Copy the volcano plot directory to the current directory
echo "Copying volcano plot to the current directory"
cp -r /opt/volcano_plot .

# Copy the differential expression file into the 'volcano_plot' directory
cp "$DIFF_EXPR_FILE" ./volcano_plot/

# Change directory into 'volcano_plot'
cd volcano_plot

# Create a unique output directory based on the output basename
OUTPUT_BASENAME="${OUTPUT_FILENAME%.html}"
OUTPUT_DIR="$OUTPUT_BASENAME"

# Run the setVars.sh script from within 'volcano_plot'
echo "Running setVars.sh with arguments: '$(basename "$DIFF_EXPR_FILE")' '$X_AXIS_COLUMN' '$Y_AXIS_COLUMN' '$LABEL_COLUMN' 'chart' 'sidebar' 'volcano' '$OUTPUT_DIR/html_data' 'OUTPUT_BASENAME"

./setVars.sh "$(basename "$DIFF_EXPR_FILE")" "$X_AXIS_COLUMN" "$Y_AXIS_COLUMN" "$LABEL_COLUMN" "chart" "sidebar" "volcano" "$OUTPUT_DIR/html_data" "$OUTPUT_BASENAME" || { echo "setVars.sh failed"; exit 1; }

# Verify the presence of the generated index.html
if [[ ! -f "volcano_plot/html_data/index.html" ]]; then
    echo "index.html not found in 'volcano_plot/html_data/'."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
# Copy 'html_data' to the output directory
cp -r "volcano_plot/html_data" "$OUTPUT_DIR/"

# Rename index.html if necessary
if [[ "$OUTPUT_FILENAME" != "index.html" ]]; then
    echo "Renaming '$OUTPUT_DIR/html_data/index.html' to '$OUTPUT_DIR/html_data/$OUTPUT_FILENAME'"
    mv "$OUTPUT_DIR/html_data/index.html" "$OUTPUT_DIR/html_data/$OUTPUT_FILENAME" || { echo "Failed to rename 'index.html' to '$OUTPUT_FILENAME'."; exit 1; }
else
    echo "Output filename is 'index.html'; no renaming needed."
fi

echo "Volcano plot generated: $OUTPUT_DIR/html_data/$OUTPUT_FILENAME"

# Additional debugging information
echo "=== Directory Structure ==="
echo "Contents of 'volcano_plot':"
ls -l

echo "Contents of '$OUTPUT_DIR':"
ls -l "$OUTPUT_DIR"

echo "Contents of '$OUTPUT_DIR/html_data':"
ls -l "$OUTPUT_DIR/html_data"
echo "=========================="

# Exit successfully
exit 0