#!/usr/bin/Rscript

# Author: Fengyuan Hu
# ORFscore calculation

# REF - https://rstudio-pubs-static.s3.amazonaws.com/164330_bee325f8e8624c18bacf64603c87de7c.html
# CMD: parallel --dry-run Rscript ~/code/github/orf-discovery/script/ORFScore.R . {} ::: `ls *.txt`
# Output header: orf.name, orfscore.raw, orfscore.log2, orfscore.p.val, orf.len.truncated, orf.counts, orf.cov

#------------------FUNC------------------------
calcORFScore <- function(frameVec) {
  # Similar to Pearson's Chi-squared Goodness of Fit Test for Count Data
  # http://stattrek.com/chi-square-test/goodness-of-fit.aspx?Tutorial=AP

  # http://www.stat.wmich.edu/s216/book/node114.html
  # The chi-square test statistic  is an overall measure of how close the observed frequencies are to the expected frequencies. 

  #p <- 1/3
  #orfscore <- log2(sum((frameVec[1]-sum(frameVec)*p)^2/(sum(frameVec)*p), (frameVec[2]-sum(frameVec)*p)^2/(sum(frameVec)*p), (frameVec[3]-sum(frameVec)*p)^2/(sum(frameVec)*p)) + 1)

  if (sum(frameVec)==0) {
    orfscore <- NA
    test.pval <- NA
  } else {
    test <- chisq.test(c(frameVec[1], frameVec[2], frameVec[3]), p=c(1/3,1/3,1/3))
    orfscore <- log2(test$statistic + 1)
    test.pval <- test$p.value

    # Apply sign based on dominant frame
    if (frameVec[1] < frameVec[2] || frameVec[1] < frameVec[3]) {
      orfscore <- orfscore * -1
    }
  }
  
  return(c(orfscore, test.pval))
}

calcCoverage <- function(countsVec) {
  pos.1.vec <- seq(from = 1, to = length(countsVec), by = 3)
  cov <- sum(countsVec[pos.1.vec]>0)/length(pos.1.vec)

  return(cov)
}

doMasking <- function(countsVec) { # filter out putative peaks, mask the position comprise more than 70% reads, mask means to asign 0 value???
  r <- max(countsVec)/sum(countsVec)
  if (!is.nan(r)) {
    if (max(countsVec)/sum(countsVec) > 0.7) {
      countsVec[which(countsVec==max(countsVec))] <- 0
    }
  }

  return(countsVec)
}
#--------------END OF FUNC----------------------

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 3) {
  orf.counts.filepath=args[1]
  outpath=args[2]
  prefix=args[3]
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

  if (length(orf.countsVec.orig) %%3 != 0) {
    message("ORF ", orf.name, " not have 3*n codons?")
    output.vec <- c(orf.name, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  } else {

    orf.countsVec.orig.pos1 <- orf.countsVec.orig[seq(1, length(orf.countsVec.orig), 3)]
    orf.countsVec.orig.pos1 <- orf.countsVec.orig.pos1[-length(orf.countsVec.orig.pos1)] # remove stop codon
    pos1.vec <- c(orf.countsVec.orig.pos1[1], orf.countsVec.orig.pos1[2], orf.countsVec.orig.pos1[3], mean(c(orf.countsVec.orig.pos1[1], orf.countsVec.orig.pos1[2])), mean(orf.countsVec.orig.pos1[3:length(orf.countsVec.orig.pos1)]), median(orf.countsVec.orig.pos1[3:length(orf.countsVec.orig.pos1)]), mean(c(orf.countsVec.orig.pos1[1], orf.countsVec.orig.pos1[2]))/mean(orf.countsVec.orig.pos1[3:length(orf.countsVec.orig.pos1)]))

    orf.cov <- calcCoverage(orf.countsVec.orig[1:(length(orf.countsVec.orig)-3)]) # remove stop codon
  
    # start codons are hyper-phased; stop codons can have differnt
    # phasing or even be de-phased depending on experimental protocol
    # so, we'll ignore the first (or more, e.g. 5) coding codon 
    #(the start codon), and the last (or more e.g. 5) coding codon. 
    codon.buffer.N <- 2
    codon.buffer.nts <- codon.buffer.N*3
  
    if (length(orf.countsVec.orig) >= 15) { # at least 4 coding codons -  start codon + 3 coding codons
      orf.countsVec <- orf.countsVec.orig[(codon.buffer.nts+1):(length(orf.countsVec.orig)-codon.buffer.nts)]
    
      orf.len.effective <- length(orf.countsVec)
    
      orf.countsVec <- doMasking(orf.countsVec)
      orf.frameVec <- colSums(matrix(orf.countsVec, ncol = 3, byrow = T))
      ORFScore <- calcORFScore(orf.frameVec)

      output.vec <- c(orf.name, ORFScore, orf.len.effective, orf.cov, pos1.vec)

    } else {
      message("ORF ", orf.name, " has less than 3 coding codons?")
      output.vec <- c(orf.name, NA, NA, 0, orf.cov, pos1.vec)
    }
  }
}, mc.cores=8)

message("Writing data to file...")
ORFScore.out.df <- data.frame(matrix(unlist(ORFScore.out), ncol=12, byrow=T))
colnames(ORFScore.out.df) <- c("ORF_ID", "ORFScore", "pval", "Effective_length", "Coverage", "codon_1_1", "codon2_1", "codon3_1", "mean_codon_1_2", "mean_codon_3_rest", "median_codon_3_rest", "ratio_mean_mean")

ORFScore.out.df$qval <- p.adjust(ORFScore.out.df$pval, "fdr")
ORFScore.out.df <- ORFScore.out.df[c("ORF_ID", "ORFScore", "pval", "qval", "Effective_length", "Coverage", "codon_1_1", "codon2_1", "codon3_1", "mean_codon_1_2", "mean_codon_3_rest", "median_codon_3_rest", "ratio_mean_mean")]

write.table(ORFScore.out.df, file=paste(outpath, "/", "ORFScore_", prefix, ".tsv", sep=""), row.names=FALSE, col.names = F, quote = F, sep = "\t")
message("Job done.")

### mocking
#countsVec.tx1 <- c(10,8,3,9,2,2,0,3,5)
#countsVec.tx2 <- c(10,10,10,10,10,10,10,10,10)
#countsVec.tx3 <- c(100,0,0,90,2,2,110,3,3)
#countsVec.tx4 <- c(100,0,0,1,0,0,0,0,0,21,2,1)

#frameVec.tx1 <- colSums(matrix(countsVec.tx1, ncol = 3, byrow = T))
#frameVec.tx2 <- colSums(matrix(countsVec.tx2, ncol = 3, byrow = T))
#frameVec.tx3 <- colSums(matrix(countsVec.tx3, ncol = 3, byrow = T))
#frameVec.tx4 <- colSums(matrix(countsVec.tx4, ncol = 3, byrow = T))

#calcORFScore(frameVec.tx1)
#calcORFScore(frameVec.tx2)
#calcORFScore(frameVec.tx3)
#calcORFScore(frameVec.tx4)
