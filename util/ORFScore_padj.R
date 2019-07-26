#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 3) {
  wd=args[1]
  orfscores.filepath=args[2]
  prefix=args[3]
  message("Working directory: ", args[1])
  message("Input file: ", args[2])
  message("Input prefix: ", args[3])
} else if (length(args) == 2) {
  wd=args[1]
  orfscores.filepath=args[2]
  message("Working directory: ", args[1])
  message("Input file: ", args[2])
} else {
  stop("Input wrong?")
}

setwd(wd)

ORFScore.out.df <- read.delim(orfscores.filepath, header=F, sep="\t")
ORFScore.out.df$V4 <- p.adjust(ORFScore.out.df$V3, "fdr")
write.table(ORFScore.out.df, file=paste("ORFScore_", prefix, ".tsv", sep=""), row.names=FALSE, col.names = F, quote = F, sep = "\t")
