#!/bin/bash

# Check for two arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <root_files_list.txt> <histogram_names_and_args.txt>"
    exit 1
fi

ROOT_FILES_LIST="$1"
HISTOGRAM_LIST="$2"

# Check if files exist
if [ ! -f "$ROOT_FILES_LIST" ]; then
    echo "Root files list '$ROOT_FILES_LIST' not found!"
    exit 1
fi

if [ ! -f "$HISTOGRAM_LIST" ]; then
    echo "Histogram list '$HISTOGRAM_LIST' not found!"
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

# Loop through each line in the histogram list
while IFS= read -r line; do
    # Read fields into positional parameters
    set -- $line

    # Check if we have at least 8 fields
    if [ $# -ne 8 ]; then
        echo "Skipping malformed line: $line"
        continue
    fi

    HISTO="$1"
    ARG1="$2"
    ARG2="$3"
    ARG3="$4"
    ARG4="$5"
    ARG5="$6"
    ARG6="$7"
    ARG7="$8"

    echo "Processing histogram: $HISTO with extra args"

    root -l -b -q "multi_histograms_General.C(\"$FILE1\", \"$FILE2\", \"$FILE3\", \"$FILE4\", \"$ARG1\", \"$HISTO\", \"$ARG2\", \"$ARG3\", \"$ARG4\", $ARG5, $ARG6, $ARG7)"
#        root -l -b -q "multi_histograms_General.C(\"$FILE1\", \"$FILE2\", \"$FILE3\",\"$FILE4\",\"Quasi-Sphere\" ,\"$histo\",\"100cm\", \"Counts / event / bin\",\"Blip multiplicity\", 0,40, true)"




done < "$HISTOGRAM_LIST"

