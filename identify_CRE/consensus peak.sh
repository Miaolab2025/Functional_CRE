## consensus peak
## We kept the peaks in >5% of samples. 
## Besides, we removed redundancy of overlapping peaks (distance between centers < 500bp) by selecting a representative peak as the one with the strongest signal. 

#!/bin/bash

# Set variables
PEAK_LIST="./peakfiles.txt"             # Text file containing paths to peak files (one per line)
FAI_FILE="./data/hg19.fa.fai"         # FAI index file of the reference genome
OUTPUT="Dup_peaks.xls"                # Output file for duplicated peaks

# Check if required input files exist
if [[ ! -f "$PEAK_LIST" ]]; then
    echo "Error: peak list file '$PEAK_LIST' not found!"
    exit 1
fi

if [[ ! -f "$FAI_FILE" ]]; then
    echo "Error: fai file '$FAI_FILE' not found!"
    exit 1
fi

# Run the Python script
echo "Running duplicate peak detection..."
srun -c 24 python ./script/merge_peaks.py \
    -AP "$PEAK_LIST" \
    -CL "$FAI_FILE" \
    -O "$OUTPUT" \

awk 'BEGIN{FS=OFS="\t"}{max =$7 ;for (i = 8; i <= NF; i++) {if ($i > max) {max = $i}}}{if (NR!=1) {print $2,$3,$4,$1,$4-$3+1,".",max,0.001,0.05,10}}' Dup_peaks.xls > Dup_peaks.bed

# Done
echo "Done. Output written to $OUTPUT"
