#!/bin/bash

# ma_plot.sh
# Script to generate MA plots from differential expression data.

# Arguments:
# $1 - Path to the differential expression file (e.g., "report.tsv")
# $2 - Name of column in file for the plot's x-axis (e.g., "baseMean")
# $3 - Name of column in file for the plot's y-axis (e.g., "log2FoldChange")
# $4 - Name of column in file for each data point's 'name' (e.g., "GeneId")
# $5 - Desired output HTML filename (e.g., "_1.html")

echo "Copying volcano plot to the current directory"
cp -r /opt/volcano_plot .
cd ./volcano_plot

# Copy the data file to the current directory
cp "$1" .

# Extract the output filename from the 5th argument
OUTPUT_FILENAME="$5"

# Derive the base name without the .html extension
OUTPUT_BASENAME=$(basename "$OUTPUT_FILENAME" .html)

# Create a unique output directory based on the output filename
OUTPUT_DIR="MD-MA_plot_${OUTPUT_BASENAME}"

# Copy the template plot directory to the unique output directory
cp -r MD-MA_plot "$OUTPUT_DIR"

# Call the variable setting script with the new output directory
./MA_PLOT_set_vars.sh "$(basename "$1")" "$2" "$3" "$4" "chart" "sideForm" "MA" "$OUTPUT_DIR/html_data"

# Rename index.html to the desired output filename
mv "$OUTPUT_DIR/html_data/index.html" "$OUTPUT_DIR/html_data/$OUTPUT_FILENAME"

echo "MA plot generated: $OUTPUT_DIR/html_data/$OUTPUT_FILENAME"