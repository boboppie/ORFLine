#!/usr/bin/Rscript

# Author: Fengyuan Hu
# ORFscore calculation
# REF - https://rstudio-pubs-static.s3.amazonaws.com/164330_bee325f8e8624c18bacf64603c87de7c.html
# CMD: parallel --dry-run Rscript ~/code/github/orf-discovery/script/ORFScore.R {} ::: `ls *.txt`
# Output header: orf.name, orfscore.raw, orfscore.log2, orfscore.p.val, orf.len.truncated, orf.counts, orf.cov

#------------------FUNC------------------------
calcORFScore <- function(frameVec) {
  # Pearson's Chi-squared Test for Count Data
  test <- chisq.test(c(frameVec[1], frameVec[2], frameVec[3]), p=c(1/3,1/3,1/3))
  orfscore.p.val <- test$p.value # no need to cast to vector as.vector()
  orfscore.raw <- test$statistic
  orfscore.log2 <- log2(orfscore.raw)
  
  # Apply sign based on dominant frame
  if (frameVec[1] < frameVec[2] || frameVec[1] < frameVec[3]) {
    orfscore.raw <- orfscore.raw * -1
    orfscore.log2 <- orfscore.log2 * -1
  }
  
  orfscoreVec <- c(orfscore.raw, orfscore.log2, orfscore.p.val)
  
  return(orfscoreVec)
}
#------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 2) {
  orf.counts.filepath=args[1]
  prefix=args[2]
  message("Input file: ", args[1])
  message("Input prefix: ", args[2])
} else if (length(args) == 1) {
  orf.counts.filepath=args[1]
  message("Input: ", args)
} else {
  stop("Input wrong?")
}

#if (grepl("/",orf.counts.filepath)) {
#  orf.filename <- unlist(strsplit(orf.counts.filepath, '/'))[-1]
#} else {
#  orf.filename <- orf.counts.filepath
#}

#orf.name <- substr(orf.filename, 1, nchar(orf.filename)-4)
#message("orf.name: ", orf.name)

#df.orf.counts.orig <- read.table(orf.counts.filepath, header = F)
#orf.countsVec.orig <- df.orf.counts.orig$V1
#message("orf.countsVec: \n")
#orf.countsVec

message("Loading data...")
no_col <- max(count.fields(orf.counts.filepath, sep = "\t"))

orf.counts <- read.table(orf.counts.filepath, sep="\t", fill=TRUE, col.names=1:no_col, header=F)
message("Data loaded.")

library(parallel)

message("Calculating ORFScore...")
ORFScore.out <- mclapply(1:nrow(orf.counts), function(i) {
  orf.count.nonNa <- orf.counts[i,which(!is.na(orf.counts[i,1:no_col]))]
  orf.name <- as.vector(t(orf.count.nonNa[1]))
  orf.countsVec.orig <- as.vector(t(orf.count.nonNa[2:length(orf.count.nonNa)]))
  
  
  # start codons are hyper-phased; stop codons can have differnt
  # phasing or even be de-phased depending on experimental protocol
  # so, we'll ignore 1 (or more, e.g. 5) codons (the start codon), 
  # and 1 (or more e.g. 5) (the stop codon).
  codon.buffer.N <- 1
  codon.buffer.nts <- codon.buffer.N*3
  
  if (length(orf.countsVec.orig) >= 9) { # 3 codons long?
    orf.countsVec <- orf.countsVec.orig[(codon.buffer.nts+1):(length(orf.countsVec.orig)-codon.buffer.nts)]
    
    orf.len <- length(orf.countsVec)
    orf.counts <- sum(orf.countsVec)
    orf.cov <- length(orf.countsVec[orf.countsVec > 0])/orf.len
    #matrix(orf.countsVec, ncol = 3, byrow = T)
    
    # to filter out putative articactual peaks, mask the most abundant position that comprised more than 70% of the total reads in the ORF
    #orf.counts.ratio <- orf.countsVec/sum(orf.countsVec)
    #orf.counts.ratio
    #which(orf.counts.ratio > 0.7)
    
    #message("orf.frameVec: \n")
    #orf.frameVec
    
    if (orf.cov == 0 | length(orf.countsVec) %%3 != 0) {
      ORFScore.vec <- c(NA, NA, NA)
    } else {
      orf.frameVec <- colSums(matrix(orf.countsVec, ncol = 3, byrow = T))
      ORFScore.vec <- calcORFScore(orf.frameVec)
    }
    
    output.vec <- c(orf.name, ORFScore.vec, orf.len, orf.counts, orf.cov)
  } else {
    message("ORF ", orf.name, " has less than 3 codons?")
    output.vec <- c(orf.name, NA, NA, NA, 0, NA, NA)
  }
}, mc.cores=8)

message("Writing data to file...")
ORFScore.out.df <- data.frame(matrix(unlist(ORFScore.out), ncol=7, byrow=T))
write.table(ORFScore.out.df, file=paste("ORFScore_", prefix, ".tsv", sep=""), row.names=FALSE, col.names = F, quote = F, sep = "\t")
message("Job done.")

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
