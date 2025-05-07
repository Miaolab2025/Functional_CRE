## This script provides a complete pipeline for processing Hi-C data using **HiC-Pro** and **Cooler**, from FASTQ files to multi-resolution `.mcool` matrices for visualization.

#!/bin/bash
# conda install -c bioconda hic-pro
# conda install -c conda-forge -c bioconda cooler

# Step 1: Build the Bowtie2 index from the reference genome
bowtie2-build {reference_genome.fa} {reference_genome}

# Step 2: Run HiC-Pro for alignment, filtering, and contact map generation
HiC-Pro -i {fastq_dir} -o {output_dir} -c config-hicpro.txt


# Step 3: Convert HiC-Pro validPairs file to .pairs format (required by cooler)
awk '{print $2"\t"$3"\t"$6"\t"$7"\t"$4"\t"$8}' {output_dir}/hicpro.validPairs > {output_dir}/hic.pairs

# Step 4: Create chrom.sizes file from the reference genome
samtools faidx {reference_genome.fa}
cut -f1,2 {reference_genome.fa.fai} > {output_dir}/chrom.sizes

# Step 5: Generate a .cool file using 10kb bin size
cooler cload pairs \
  --assembly {output_dir}/chrom.sizes \
  -c1 2 -p1 3 -c2 4 -p2 5 \
  chrom.sizes:1000 {output_dir}/hic.pairs {output_dir}/hic_1kb.cool

# Step 6: Create a multi-resolution .mcool file for visualization
cooler zoomify -r 10000,100000 -o {output_dir}/hic.mcool {output_dir}/hic_1kb.cool












