#!/usr/bin/Rscript

# Author: Fengyuan Hu
# ORFscore calculation
# REF - https://rstudio-pubs-static.s3.amazonaws.com/164330_bee325f8e8624c18bacf64603c87de7c.html
# CMD: parallel --dry-run Rscript ~/code/github/orf-discovery/script/ORFScore.R {} ::: `ls *.txt`

calcORFScore <- function(frameVec) {
  orfscore.raw <- as.vector(chisq.test(c(frameVec[1], frameVec[2], frameVec[3]),
                                   p=c(1/3,1/3,1/3))$statistic)
  
  orfscore.log2 <- log2(as.vector(chisq.test(c(frameVec[1], frameVec[2], frameVec[3]),
                                        p=c(1/3,1/3,1/3))$statistic))
  
  # Apply sign based on dominant frame
  if (frameVec[1] < frameVec[2] || frameVec[1] < frameVec[3]) {
    orfscore.raw <- orfscore.raw * -1
    orfscore.log2 <- orfscore.log2 * -1
  }
  
  orfscoreVec <- c(orfscore.raw, orfscore.log2)
  
  return(orfscoreVec)
}

args <- commandArgs(trailingOnly = TRUE)
message("Input: ", args)

if(length(args) == 1) {
  orf.counts.filepath=args[1]
} else {
  stop("Input wrong?")
}

if (grepl("/",orf.counts.filepath)) {
  orf.filename <- unlist(strsplit(orf.counts.filepath, '/'))[-1]
} else {
  orf.filename <- orf.counts.filepath
}

orf.name <- substr(orf.filename, 1, nchar(orf.filename)-4)
#message("orf.name: ", orf.name)

df.orf.counts <- read.table(orf.counts.filepath, header = F)
orf.countsVec <- df.orf.counts$V1
#message("orf.countsVec: \n")
#orf.countsVec

orf.cov <- length(orf.countsVec[orf.countsVec > 0])/length(orf.countsVec)
#matrix(orf.countsVec, ncol = 3, byrow = T)

# to filter out putative articactual peaks, mask the most abundant position that comprised more than 70% of the total reads in the ORF
orf.counts.ratio <- orf.countsVec/sum(orf.countsVec)
#orf.counts.ratio
#which(orf.counts.ratio > 0.7)

#message("orf.frameVec: \n")
#orf.frameVec

if (orf.cov == 0 | length(orf.countsVec) %%3 != 0) {
  ORFScore.vec <- c(NA, NA)
} else {
  orf.frameVec <- colSums(matrix(orf.countsVec, ncol = 3, byrow = T))
  ORFScore.vec <- calcORFScore(orf.frameVec)
}

output.vec <- c(orf.name, ORFScore.vec, orf.cov)

dir.create(file.path(".", "ORFScore"), showWarnings = FALSE)
setwd(file.path(".", "ORFScore"))
write.table(t(output.vec), file =paste(orf.name, ".ORFScore.tsv", sep = ""), row.names=FALSE, col.names = F, quote = F, sep = "\t")

### mocking
#countsVec.tx1 <- c(10,8,3,9,2,2,0,3,5)
#countsVec.tx2 <- c(10,10,10,10,10,10,10,10,10)
#countsVec.tx3 <- c(100,0,0,90,2,2,110,3,3)
#countsVec.tx4 <- c(100,0,0,1,0,0,0,0,0,21,2,1)

#frameVec.tx1 <- colSums(matrix(countsVec.tx1, nrow = 3, byrow = T))
#frameVec.tx2 <- colSums(matrix(countsVec.tx2, nrow = 3, byrow = T))
#frameVec.tx3 <- colSums(matrix(countsVec.tx3, nrow = 3, byrow = T))
#frameVec.tx4 <- colSums(matrix(countsVec.tx4, nrow = 3, byrow = T))

#calcORFScore(frameVec.tx1)
#calcORFScore(frameVec.tx2)
#calcORFScore(frameVec.tx3)
#calcORFScore(frameVec.tx4)
