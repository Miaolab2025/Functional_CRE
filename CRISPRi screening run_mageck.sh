## CRISPRi screening for functional CREs with MAGeCK
## In this part, the number of reads for each sgRNA were counted and normalized to the total reads of all sgRNAs using MAGeCK identify potential functional CREs positively or negatively enriched in the screening.  

#!/bin/bash
# Step 1: Count and normalize sgRNA Reads
mageck count \
  -k sample_counts.txt \
  -t Treatment_Rep1,Treatment_Rep2 \
  -c Control_Rep1,Control_Rep2 \
  --norm-method median \
  -o output_prefix

# Step 2: Test for differential sgRNA abundance
mageck test \
  -k output_prefix.count.txt \
  -t Treatment_Rep1,Treatment_Rep2 \
  -c Control_Rep1,Control_Rep2 \
  -n output_prefix \
  --adjust-method fdr \
  --control-sgrna control_sgrnas.txt \
  --gene-lfc-method median
