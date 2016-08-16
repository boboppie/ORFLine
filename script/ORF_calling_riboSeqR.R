#!/usr/bin/Rscript

# Author: Fengyuan Hu
# Using riboSeqR to call ORFs
# Input:
# Transcriptome alignment file in BAM format
# Transcriptome (protein coding + lincRNAs) in Fasta format

suppressPackageStartupMessages(library(riboSeqR))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))

setwd("~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-18_10-44-28/ORFCalling/riboSeqR")

# Find novel CDS
# The findCDS function defines (from a fasta file of the transcriptome) all possible coding sequences with a start and stop codon in frame.
refTranscriptomeFile <- "~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/pc_and_lncRNA_transcripts.fa"
#refGenomeForwardFile <- "~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa"
# refGenomeReverseFile <- "~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.revcomp.fa"
#refGenomeFai.df <- read.delim("~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa.fai", header = F)

allORFs.transcriptome <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("ATG", "CTG", "TTG", "GTG"), stopCodon = c("TAG", "TAA", "TGA"))
#allORFs.genome.forward <- findCDS(fastaFile = refGenomeForwardFile, startCodon = c("ATG", "CTG", "TTG", "GTG"), stopCodon = c("TAG", "TAA", "TGA"))
#allORFs.genome.reverse <- findCDS(fastaFile = refGenomeReverseFile, startCodon = c("ATG", "CTG", "TTG", "GTG"), stopCodon = c("TAG", "TAA", "TGA"))

###### TESTING 
#refTranscriptomeFile <- "~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/test/Cxcr4.fa"
refTranscriptomeFile <- "~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/test/Ptbp3.fa"
allORFs.transcriptome <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
######

# remove the bit after space in seqnames, e.g. "chr1 1" -> "chr1"
#seqlevels(allORFs.genome.forward) <- gsub(" .*$", "",seqlevels(allORFs.genome.forward))
#seqlevels(allORFs.genome.reverse) <- gsub(" .*$", "",seqlevels(allORFs.genome.reverse))

#strand(allORFs.genome.forward) <- rep("+", length(allORFs.genome.forward))
#strand(allORFs.genome.reverse) <- rep("-", length(allORFs.genome.reverse))

#seqlengths(allORFs.genome.forward) <- refGenomeFai.df$V2
#seqlengths(allORFs.genome.reverse) <- refGenomeFai.df$V2

#save(allORFs.genome.forward, file="allORFsGenome.forward.rda")
#save(allORFs.genome.reverse, file="allORFsGenome.reverse.rda")

#allORFs.genome.reverse.orig <- allORFs.genome.reverse
# TODO create new IRange and replace the old one?
#CORES <- 4
#newranges <- mclapply(1:length(allORFs.genome.reverse), function(i) {
#  chr.len <- seqlengths(allORFs.genome.reverse[i])[runValue(seqnames(allORFs.genome.reverse[i]))][[1]]
#  newstart = chr.len - end(allORFs.genome.reverse[i]) + 1
#  newend = chr.len - start(allORFs.genome.reverse[i]) + 1
#  c(newstart, newend)
#}, mc.cores=CORES)

#newranges.df <- as.data.frame(do.call(rbind, newranges))
#newirange <- IRanges(start = newranges.df$V1, end = newranges.df$V2)
#allORFs.genome.reverse@ranges <- newirange

#save(allORFs.genome.reverse, file="allORFsGenome.reverse.new.rda")

#allORFs.genome.allFrames.gr <- c(allORFs.genome.forward, allORFs.genome.reverse)

#save(allORFs.genome.allFrames.gr, file="allORFsGenome.rda")
#load("allORFsGenome.rda")
#export(allORFs.genome.allFrames.gr, "allORFsGenome.bed")

# Then count the reads and their frame mappings
riboDatFile.uniqHit <- "~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-18_10-44-28/bowtie-transcriptome/merged_q255.bam"
riboDatFile.multiHit <- "~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-18_10-44-28/bowtie-transcriptome/merged_multiple.bam"

riboDat.uniqHit <- readRibodata(riboDatFile.uniqHit, replicates = "wt1", zeroIndexed = F)
riboDat.multiHit <- readRibodata(riboDatFile.multiHit, replicates = "wt1", zeroIndexed = F)

len <- c(25:35)
fCs.allORFs.transcriptome.uniqHit <- frameCounting(riboDat.uniqHit, allORFs.transcriptome, lengths = len)
fS.allORFs.transcriptome.uniqHit <- readingFrame(rC = fCs.allORFs.transcriptome.uniqHit, lengths = len); fS.allORFs.transcriptome.uniqHit
plotFS(fS.allORFs.transcriptome.uniqHit, lengths = len)

fCs.allORFs.transcriptome.multiHit <- frameCounting(riboDat.multiHit, allORFs.transcriptome, lengths = len)
fS.allORFs.transcriptome.multiHit <- readingFrame(rC = fCs.allORFs.transcriptome.multiHit, lengths = len); fS.allORFs.transcriptome.multiHit
plotFS(fS.allORFs.transcriptome.multiHit, lengths = len)

len.RPF <- c(26:34)
# and filter on the number of ribosomal fragments mapping; in this data, I was expecting ribosomal fragments of length 27 and to map to frame 1 relative to the start codon
ffCs.allORFs.transcriptome.uniqHit <- filterHits(fCs.allORFs.transcriptome.uniqHit, lengths = len.RPF, frame = 0, hitMean = 50, unqhitMean = 10, fS  = fS.allORFs.transcriptome.uniqHit)


save(ffCs.allORFs.genome, file="ORFCalledFrame0.rda")

# save GRange to BED format (rtracklayer)
export(ffCs.allORFs.genome@CDS, "ORFCalledFrame0.bed")
