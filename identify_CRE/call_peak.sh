#!/bin/bash

###call peak
input_dir="./bam_files"   
output_base="./macs2_peaks"  

genome_size="hs"           
shift_size=100             
ext_size=200               
buffer_size=100000        
qval_thresh=0.01           
fe_cutoff=1                

# parameters for ChIP-seq
#genome_size="hs"             # Genome size; use "hs" for human, "mm" for mouse, or specific size like 2.7e9
#shift_size=100               # Shift size (for calculating --slocal and --llocal)
#ext_size=147                 # Fragment extension size
#mfold_min=10                 # Minimum fold-enrichment for model building
#mfold_max=50                 # Maximum fold-enrichment for model building
#buffer_size=100000           # Buffer size to reduce memory usage
#qval_thresh=0.01             # FDR q-value cutoff for peak detection
#broad_cutoff=0.01            # q-value cutoff for broad peak detection
#fe_cutoff=1                  # Fold enrichment cutoff
#tsize=150                    # Tag size (read length)
#bw=147                       # Bandwidth for fragment size estimation
#keep_dup=1                   # Number of duplicate reads to keep

mkdir -p "$output_base"

for bam_file in "$input_dir"/*.bam; do
    filename=$(basename "$bam_file")
    sample_name="${filename%%.*}"
    sample_outdir="${output_base}/${sample_name}"
    mkdir -p "$sample_outdir"
    echo "Processing $sample_name..."
    macs2 callpeak \
        -t "$bam_file" \
        -g "$genome_size" \
        -n "$sample_name" \
        --outdir "$sample_outdir" \
        --buffer-size "$buffer_size" \
        -f BAMPE \
        --shift "$shift_size" \
        --extsize "$ext_size" \
        --nomodel \
        --keep-dup all \
        --SPMR \
        -q "$qval_thresh" \
        --fe-cutoff "$fe_cutoff"
done
echo "MACS2 peak calling finished."
