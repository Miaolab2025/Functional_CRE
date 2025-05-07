## Reads mapped to the exon regions of each gene were counted by featureCounts (Subread-1.5.1; Bioconductor) for RNA-seq data.

featureCounts -T 10 \
    -a ~/orign.data/gene.annotation.file/gencode.v47lift37.annotation.gtf.gz \
    -o gene_counts.txt \
    -p -B -C -g gene_id \
    *.bam

head -n 1 gene_counts.txt | cut -f7- > header.txt
tail -n +2 gene_counts.txt | cut -f1,7- > counts_body.txt
paste <(echo -e "Geneid\n$(cut -f1 counts_body.txt)") counts_body.txt | cut -f1,2- > read_count_matrix.txt

