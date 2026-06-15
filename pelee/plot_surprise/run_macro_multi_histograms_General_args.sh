#!/bin/bash

# Check for one or two arguments
#if [ $# -lt 1 ] || [ $# -gt 2 ]; then
#    echo "Usage: $0 <root_files_list.txt> [plot_configs.txt]"
#    exit 1
#fi

ROOT_FILES_LIST="${1:-files.txt}"   # default to files.txt if not provided
PLOT_CONFIG="${2:-Splot_configs.txt}"   # default to plot_configs.txt

# Check if files exist
if [ ! -f "$ROOT_FILES_LIST" ]; then
    echo "Root files list '$ROOT_FILES_LIST' not found!"
    exit 1
fi

if [ ! -f "$PLOT_CONFIG" ]; then
    echo "Plot configuration file '$PLOT_CONFIG' not found!"
    exit 1
fi

# Read the 4 ROOT file paths
count=0
while IFS= read -r line; do
    case $count in
        0) FILE1="$line" ;;
        1) FILE2="$line" ;;
        2) FILE3="$line" ;;
        3) FILE4="$line" ;;
        *) echo "Error: More than 4 root files found in $ROOT_FILES_LIST" >&2; exit 1 ;;
    esac
    count=$((count + 1))
done < "$ROOT_FILES_LIST"

if [ $count -ne 4 ]; then
    echo "Error: Expected exactly 4 root files, but found $count"
    exit 1
fi

echo "Processing all histogram variable sets..."
echo "Input files:"
echo "  File 1: $FILE1"
echo "  File 2: $FILE2" 
echo "  File 3: $FILE3"
echo "  File 4: $FILE4"
echo "  Config: $PLOT_CONFIG"

# Run the ROOT macro once - it will process all variable sets from config file
root -l -b -q "surprise_multi_histograms_General.C(\"$FILE1\", \"$FILE2\", \"$FILE3\", \"$FILE4\", \"$PLOT_CONFIG\")" 

echo "All plots saved to plots_output"
