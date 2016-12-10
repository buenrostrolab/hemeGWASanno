
#+ echo=FALSE, message=FALSE, warning=FALSE
library(readr)
library(plotly)
library(heatmaply)
library(scTools)
library(BuenColors)
library(ggplot2)
library(data.table)


if (basename(getwd()) != "code") setwd("code")

# hemecounts <- data.matrix(fread("../data/Corces/GSE74912_ATACseq_All_Counts.txt"))
# h <- hemecounts[,4:80]
# cc <- colnames(h)
# shortnames <- sapply(1:length(cc), function(s){
#   if(s <= 23 | (s <= 74 & s >= 62 )){ strsplit(cc[s], split = "-")[[1]][2]
#   } else { strsplit(cc[s], split = "-")[[1]][3]}
# })
# noa <- gsub("A", "", shortnames)
# nob <- gsub("B", "", noa)
# noc <- gsub("C", "", nob)
# 
# n <- as.character(seq(1:15))
# cell <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "missing", "CD4", "CD8", "NK", "missing", "Bcell", "CLP", "Ery")
# names(cell) <- n
# colnames(h) <- paste0(unname(cell[noc]), "_",1:length(noc))
# 
# type <- read.table("../data/Corces/type.atacpeak.txt")
# cor.TSS <- cor(h[type[,1] =="promoter-TSS",], use="complete.obs", method="pearson")
# saveRDS(cor.TSS, file = "../data/Corces/cor.TSS.rds")
# 
# cor.outs <- cor(h[type[,1] =="Intergenic",], use="complete.obs", method="pearson")
# saveRDS(cor.outs, file = "../data/Corces/cor.outs.rds")
# 
# 
cor.outs <- readRDS("../data/Corces/cor.outs.rds")
cor.TSS <- readRDS("../data/Corces/cor.TSS.rds")

#' # All versus All
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
mco <- melt(cor.outs)
mct <- melt(cor.TSS)

plotdf <- data.frame(mco, mct)
plotdf <- plotdf[plotdf$value < 1, c(1,2,3,6)]
names(plotdf) <- c("Name1", "Name2", "Distal", "Promoter")
plotdf$anno <- apply(sapply(colnames(plotdf), function(name){paste0(name, ": ", plotdf[,name])}), 1, paste, collapse = "<br>")
plotdf$Difference <- plotdf$Promoter - plotdf$Distal

names(plotdf) <- c("Name1", "Name2", "Distal", "Promoter", "anno", "Difference")
plot_ly(plotdf, x = ~Promoter, y = ~Distal, mode = "markers", color = ~Difference, text = ~anno,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) 

# Slow loop
p1 <- c("HSC", "MPP", "MPP", "LMPP","LMPP","CLP", "CLP", "CLP", "CLP", "CMP", "CMP", "MEP", "GMP")
p2 <- c("MPP","LMPP", "CMP", "CLP", "GMP", "CD4", "CD8", "Bcell",   "NK", "GMP",  "MEP", "Ery", "Mono")
diffpart <- apply(plotdf, 1, function(t){
  ct1 <- strsplit(t[1], split = "_")[[1]][1]
  ct2 <- strsplit(t[2], split = "_")[[1]][1]
  boo <- 0
  for(ii in 1:length(p1)){
    if(p1[ii] == ct1 & p2[ii] == ct2) boo <- 1
    if(p1[ii] == ct2 & p2[ii] == ct1) boo <- 2
  }
  boo
})

plotdf$parent <- diffpart
plotdf2 <- plotdf[unname(diffpart) == 1 | unname(diffpart) == 2, ]
plotdf2$child <- plotdf2$parent*(-1) + 3

#' # Only pairs that are related
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
plot_ly(plotdf2, x = ~Promoter, y = ~Distal, mode = "markers", color = ~Difference, text = ~anno,
         colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), marker = list(size = 10)) 

#' # Colored by parent
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
parentCell <- sapply(1:dim(plotdf2)[1], function(rowidx){
  strsplit(as.character(plotdf2[rowidx,plotdf2[rowidx,"parent"]]),split = "_")[[1]][1]
})
plotdf2$parentCell <- parentCell

plot_ly(plotdf2, x = ~Promoter, y = ~Distal, mode = "markers", color = ~parentCell, text = ~anno,
         marker = list(size = 10)) 

#' # Colored by parent:child
#+ echo=FALSE, message=FALSE, warning=FALSE, fig.height = 7, fig.width = 7
parentChild <- sapply(1:dim(plotdf2)[1], function(rowidx){
  paste0(strsplit(as.character(plotdf2[rowidx,plotdf2[rowidx,"parent"]]),split = "_")[[1]][1],
         ":", strsplit(as.character(plotdf2[rowidx,plotdf2[rowidx,"child"]]),split = "_")[[1]][1])
})

plot_ly(plotdf2, x = ~Promoter, y = ~Distal, mode = "markers", color = ~parentChild, text = ~anno,
         marker = list(size = 10)) 