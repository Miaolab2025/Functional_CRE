## --------------------------
## Identify differential CREs during pairwise comparison
## Data: ATAC-seq or H3K27ac ChIP-seq  (example: ATAC-seq)
## Comparisons: User specifies input files 
## --------------------------

library(csaw)
library(edgeR)
library(rtracklayer)
library(ChIPseeker)

# Step 1: Input BAM files
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide exactly two input BAM file directories for comparison.")
}

bam.dir1 <- args[1]  # First comparison group 
bam.dir2 <- args[2]  # Second comparison group

# Read BAM files for each group
bam.files1 <- list.files(bam.dir1, pattern = "\\.bam$", full.names = TRUE)
bam.files2 <- list.files(bam.dir2, pattern = "\\.bam$", full.names = TRUE)

# Flatten BAM list
bam.list <- c(bam.files1, bam.files2)
group.labels <- rep(c("Group1", "Group2"), times = c(length(bam.files1), length(bam.files2)))

# Step 2: Preprocessing of input file
param <- readParam(
  dedup = TRUE,
  minq = 20
)

# Sliding Window Counting
window.width <- 150
spacing <- 50
win.data <- windowCounts(
  bam.list,
  ext = 200,
  spacing = spacing,
  width = window.width,
  param = param
)

# Filter low-abundance windows
keep <- aveLogCPM(asDGEList(win.data)) > 0.5
win.data <- win.data[keep,]

# Normalization
win.data <- normFactors(win.data, se.out = TRUE)

# Step 3: Design matrix for two-group comparison
design <- model.matrix(~0 + factor(group.labels))
colnames(design) <- c("Group1", "Group2")

# ---------- edgeR model fitting ----------
y <- asDGEList(win.data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

# ---------- Define contrast ----------
contrast <- makeContrasts(Group2 - Group1, levels = design)

# ---------- Run test and merge peaks ----------
res <- glmQLFTest(fit, contrast = contrast)
merged <- mergeWindows(rowRanges(win.data), tol = 2000)
combined <- combineTests(merged$id, res$table)

sig <- combined[combined$FDR < 0.05 & abs(combined$logFC) > 1, ]

# Export to BED
bed.file <- paste0("diff_peaks_", basename(bam.dir1), "_vs_", basename(bam.dir2), ".bed")
export.bed(sig, con = bed.file)

# Optional: Peak annotation
# annotation <- annotatePeak("diff_peaks.bed", tssRegion = c(-3000, 3000), TxDb = txdb)
# plotAnnoPie(annotation)
