## Characterization of differential CREs
## We conducted positional annotation  on the differential CREs using ChIPseeker

install.packages('BiocManager')
BiocManager::install('ChIPseeker')
BiocManager::install("GenomicFeatures")
BiocManager::install("org.Hs.eg.db")


library(ChIPseeker)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Read the differential peak file (BED format)
peak.file <- "diff_peaks.bed"  # ← Replace with your own BED file
peaks <- readPeakFile(peak.file)

# Perform peak annotation (e.g., within ±3kb of TSS, using hg19 genome)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(peaks,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")

# Save annotation result to a tab-delimited file
write.table(as.data.frame(peakAnno), 
            file = sub(".bed$", ".annotated.txt", peak.file), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


