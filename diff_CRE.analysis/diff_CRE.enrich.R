## In this part we characterized the CREs through enrichment analyses,

# Load required packages
library(data.table)
library(ggplot2)
library(devtools)
install_github("wjawaid/enrichR")  # Install the latest enrichR package from GitHub
library(enrichR)
library(readxl)

# Set enrichR site for human
setEnrichrSite("Enrichr")

# List all available Enrichr databases
dbs_available <- listEnrichrDbs()

# Select databases to be used
dbs <- c(
  "GO_Biological_Process_2024",
  "GO_Cellular_Component_2024",
  "GO_Molecular_Function_2024",
  "KEGG_2021_Human"
)

# Set working directory
setwd("~/path/to/diff_CRE_annotation")

# Load annotated differential peak file
peakAnno <- fread("diff_peaks.annotated.txt", header = TRUE)

# Extract gene symbols from CREs
genes_CRE     <- peakAnno[type == "CRE" & !is.na(SYMBOL), SYMBOL]

# Perform enrichment analysis for CRE_gene lists
eCRE    <- enrichr(genes_CRE, dbs)

# Extract KEGG enrichment results
CRE_kegg    <- eCRE[["KEGG_2021_Human"]]

# Add group labels
CRE_kegg$type    <- "CRE"

# Select top 20 enriched pathways by Combined Score
CRE_kegg    <- CRE_kegg[1:20, ][order(Combined.Score), ]

# Convert Term column to factor to control plot order
CRE_kegg$Term <- factor(CRE_kegg$Term, levels = CRE_kegg$Term)

# Export KEGG enrichment results
write.table(CRE_kegg, file = "KEGG_enrichment_CRE.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
