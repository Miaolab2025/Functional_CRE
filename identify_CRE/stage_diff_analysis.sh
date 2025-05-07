## In this part, to identify significantly differential cis-regulatory elements (CREs) during the transition from normal tissue to advanced adenoma and cancer, we performed pairwise differential peak analysis of H3K27ac ChIP-seq or ATAC-seq data using the csaw package. 
## Each comparison (C vs. N, C vs. A, and A vs. N) was conducted independently by specifying two input sample groups in the stage_diff_analysis.sh script.

#!/bin/bash

# Define the directories for the two input groups   (example: ATAC-seq, C vs N)
GROUP1_DIR="/path/to/C"
GROUP2_DIR="/path/to/N"

# Run the R script with the input directories as arguments
Rscript differential_CRE_analysis.R $GROUP1_DIR $GROUP2_DIR
