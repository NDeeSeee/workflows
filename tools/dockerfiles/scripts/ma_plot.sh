#!/bin/bash

# ma_plot.sh
# Script to generate MA plots from differential expression data.

set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error.
set -o pipefail  # Catch errors in pipelines.

#########################
# Function Definitions  #
#########################

# Function to print error messages and exit
error_exit() {
    echo "Error: $1" >&2
    exit 1
}

# Function to print debugging information
print_debug_info() {
    echo "=== ma_plot.sh Debugging Information ==="
    echo "Argument 1 (diff_expr_file.path): $DIFF_EXPR_FILE"
    echo "Argument 2 (x_axis_column): $X_AXIS_COLUMN"
    echo "Argument 3 (y_axis_column): $Y_AXIS_COLUMN"
    echo "Argument 4 (label_column): $LABEL_COLUMN"
    echo "Argument 5 (output_filename): $OUTPUT_FILENAME"
    echo "========================================="
}

# Function to copy files and directories with appropriate flags
copy_file_or_dir() {
    local src="$1"
    local dest="$2"
    echo "Copying '$src' to '$dest'"

    if [[ -d "$src" ]]; then
        cp -r "$src" "$dest" || error_exit "Failed to copy directory '$src' to '$dest'."
    elif [[ -f "$src" ]]; then
        cp "$src" "$dest" || error_exit "Failed to copy file '$src' to '$dest'."
    else
        error_exit "Source '$src' is neither a file nor a directory."
    fi
}

# Function to check file existence
check_file_exists() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        error_exit "File '$file' does not exist."
    fi
}

#########################
# Argument Parsing      #
#########################

# Ensure exactly 5 arguments are provided
if [[ $# -ne 5 ]]; then
    error_exit "Expected 5 arguments, but got $#.
Usage: ma_plot.sh <diff_expr_file> <x_axis_column> <y_axis_column> <label_column> <output_filename>"
fi

# Assign arguments to variables
DIFF_EXPR_FILE="$1"
X_AXIS_COLUMN="$2"
Y_AXIS_COLUMN="$3"
LABEL_COLUMN="$4"
OUTPUT_FILENAME="$5"

# Print debugging information
print_debug_info

#########################
# Main Script Execution #
#########################

# Define constants
VOLCANO_PLOT_SOURCE="/opt/volcano_plot"
VOLCANO_PLOT_DEST="volcano_plot"

# Copy volcano plot directory
copy_file_or_dir "$VOLCANO_PLOT_SOURCE" "$VOLCANO_PLOT_DEST"

# Change to the volcano_plot directory
cd "$VOLCANO_PLOT_DEST" || error_exit "Failed to change directory to '$VOLCANO_PLOT_DEST'."

# Copy the differential expression file
copy_file_or_dir "$DIFF_EXPR_FILE" "./"

# Extract the basename of the differential expression file
DATA_FILE="$(basename "$DIFF_EXPR_FILE")"

# Verify the copied differential expression file exists
check_file_exists "$DATA_FILE"

echo "Successfully copied '$DATA_FILE' to the current directory."

# Derive the base name without the .html extension
OUTPUT_BASENAME="${OUTPUT_FILENAME%.html}"

echo "Output basename: $OUTPUT_BASENAME"

# Create a unique output directory based on the output basename
OUTPUT_DIR="MD-MA_plot_${OUTPUT_BASENAME}"

echo "Creating output directory: $OUTPUT_DIR"
copy_file_or_dir "MD-MA_plot" "$OUTPUT_DIR" || error_exit "Failed to copy MD-MA_plot to '$OUTPUT_DIR'."

# Verify the presence of the variable setting script
SET_VARS_SCRIPT="./MA_PLOT_set_vars.sh"
if [[ ! -x "$SET_VARS_SCRIPT" ]]; then
    error_exit "MA_PLOT_set_vars.sh script not found or not executable."
fi

# Call the variable setting script with the required arguments
echo "Running variable setting script: $SET_VARS_SCRIPT"
"$SET_VARS_SCRIPT" "$DATA_FILE" "$X_AXIS_COLUMN" "$Y_AXIS_COLUMN" "$LABEL_COLUMN" "chart" "sideForm" "MA" "$OUTPUT_DIR/html_data" || error_exit "MA_PLOT_set_vars.sh failed."

# Verify the presence of the generated index.html
GENERATED_INDEX="$OUTPUT_DIR/html_data/index.html"
if [[ ! -f "$GENERATED_INDEX" ]]; then
    error_exit "index.html not found in '$OUTPUT_DIR/html_data/'."
fi

# Determine if renaming is needed
if [[ "$OUTPUT_FILENAME" != "index.html" ]]; then
    echo "Renaming '$GENERATED_INDEX' to '$OUTPUT_DIR/html_data/$OUTPUT_FILENAME'"
    mv "$GENERATED_INDEX" "$OUTPUT_DIR/html_data/$OUTPUT_FILENAME" || error_exit "Failed to rename '$GENERATED_INDEX' to '$OUTPUT_FILENAME'."
else
    echo "Output filename is 'index.html'; no renaming needed."
fi

echo "MA plot generated: $OUTPUT_DIR/html_data/$OUTPUT_FILENAME"

# Additional debugging information
echo "=== Directory Structure ==="
echo "Current directory contents:"
ls -l

echo "Contents of '$OUTPUT_DIR':"
ls -l "$OUTPUT_DIR"

echo "Contents of '$OUTPUT_DIR/html_data':"
ls -l "$OUTPUT_DIR/html_data"
echo "=========================="

# Exit successfully
exit 0