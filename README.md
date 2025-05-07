# Functional_CRE
This repository contains the code used to perform the analysis in the "Combined epigenomics and CRISPRi screening characterize dynamics of cis-regulatory elements and functional variants to promote colorectal cancer development" paper. The purpose of making the code available is for transparency and data analysis reproducibility.

# Overview
Genetic variants associated with colorectal cancer (CRC) are predominantly noncoding and often located in regions of high linkage disequilibrium, particularly within cis-regulatory elements (CREs). Dissecting the causal variants and uncovering their underlying molecular mechanisms remains a major challenge. To address this, Functional_CRE offers a modular and integrative pipeline for the identification, annotation, and functional interpretation of CREs involved in CRC development. By leveraging multi-omics data from ATAC-seq, H3K27ac ChIP-seq, Hi-C, and RNA-seq, this pipeline enables the detection of differential CREs across disease progression stages (from normal tissue to adenoma to carcinoma), the assignment of putative target genes, and the inference of distal regulatory interactions through 3D chromatin architecture. Additionally, it supports functional annotation through enrichment analyses, and incorporates CRISPRi screening data to prioritize regulatory variants with potential roles in tumorigenesis.

# Workflow
The analyses contained in our work:
	.[RNA-seq reads count featureCount.sh] - This work flow processes RNA-seq data to quantify read counts specifically over annotated exonic regions.
	.[identify_CRE] - This workflow is used to identify differential CREs, including peak calling, consensus peak identification, and stage-specific differential peak analysis, based on ATAC-seq and H3K27ac ChIP-seq data.
	.[diff_CRE.analysis] - This workflow is designed for the positional annotation and functional enrichment analysis of cis-regulatory elements (CREs), using dedicated scripts including diff_CRE.annot.R and diff_CRE.enrich.R.
	.[CRISPRi screening run_mageck.sh] - This workflow enables the systematic processing and analysis of CRISPRi screening data from cell lines, supporting the identification of functional cis-regulatory elements (CREs) involved in gene regulation.
	.[HiC_analysis] - This workflow processes Hi-C data from raw reads to contact matrices and identifies TAD boundaries based on insulation score analysis.
# License
miaolab
