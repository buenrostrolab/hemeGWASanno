#+ echo=FALSE, message=FALSE, warning=FALSE
library(chromVAR)
library(Matrix)
library(GenomicRanges)
library(diffloop)

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

load("../output/peaksCounts.rda")
types <- readRDS("../output/sampleTypes.rds")
sample_annotation <- DataFrame(type = types)
hemedf <- SummarizedExperiment(assays = list(counts = Matrix(counts)), rowRanges = peaks, colData = sample_annotation)

# Read GWAS hits in and process ||| BED FILE
hitdf <- read.table("../data/GWAS/20160617_both_sig_MPRA_SNPs.bed")
names(hitdf) <- c("chr", "start", "end", "snp")
hitG <- GRanges(hitdf)

# Read two lines
hitdf <- read.table("../data/GWAS/sarc.ea.10e-04.txt")
hitdf$V3 <- hitdf$V2
names(hitdf) <- c("chr", "start", "end")
hitG <- addchr(GRanges(hitdf))

# Build matches index
ov <- findOverlaps(peaks, hitG)
ix <- Matrix(as.integer(1:length(peaks) %in% unique(queryHits(ov))))
ixx <- SummarizedExperiment(assays = list(match = Matrix(ix)), rowRanges = peaks,  colData = DataFrame(hit = "GWAS"))

deviations <- compute_deviations(hemedf, background_peaks = bg, annotations = ixx)

rawdf <- data.frame( types = types, zscore = unname(assays(deviations)$z)[1,])
aggregate(rawdf, by = list(Cell = rawdf$types), FUN = mean)[,c(1,3)]
