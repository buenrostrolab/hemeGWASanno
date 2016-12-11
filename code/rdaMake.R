
#+ echo=FALSE, message=FALSE, warning=FALSE
library(GenomicRanges)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVAR)

if (basename(getwd()) != "code") setwd("code")

#Process counts

hemecounts <- data.matrix(fread(input = 'zcat < ../data/heme/GSE74912_ATACseq_All_Counts.txt.gz'))
h <- hemecounts[,4:80]
cc <- colnames(h)
shortnames <- sapply(1:length(cc), function(s){
   if(s <= 23 | (s <= 74 & s >= 62 )){ strsplit(cc[s], split = "-")[[1]][2]
   } else { strsplit(cc[s], split = "-")[[1]][3]}
})
noa <- gsub("A", "", shortnames)
nob <- gsub("B", "", noa)
noc <- gsub("C", "", nob)
 
n <- as.character(seq(1:15))
cell <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "missing", "CD4", "CD8", "NK", "missing", "Bcell", "CLP", "Ery")
names(cell) <- n
colnames(h) <- paste0(unname(cell[noc]), "_",1:length(noc))

#Process peaks
peakdf <- read.table("../data/heme/peaks.bed", header = TRUE)
names(peakdf) <- c("chr", "start", "end")
peaks <- GRanges(peakdf)

peaks <- get_gc(peaks, genome = BSgenome.Hsapiens.UCSC.hg19)
praw <- peaks
peaks <- sort(peaks)
counta <- h[match(peaks,praw),]
counts <- counta[rowSums(counta) != 0,]
peaks <- peaks[rowSums(counta) != 0]

#' ## Get motifs; compute deviations. 
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = FALSE
bg <- get_background_peaks(counts_mat = counts, bias = peaks)

save(counts, peaks, bg, file = "../output/peaksCounts.rda")