#!/usr/bin/Rscript

# Author: Fengyuan Hu
# Search ORFs from a given transcript sequences, covert local coordinates to genomic coordinates in BED format
# Input:
# Transcriptome in Fasta format
# Transcriptome annotation in GTF format

# Make sure annotation is consistent with sequences (same version), or extract sequences from genome?

org <- "Mus musculus"
startCodon <- "ATG"
refTranscriptomeFile <- "../../out/ref/pipeline_global/transcriptome.fa"
gtfFile <- "../../out/ref/pipeline_global/annotation.gtf"
NCORE <- 1

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 5) {
  org=args[1]
  startCodon=args[2]
  refTranscriptomeFile=args[3]
  gtfFile=args[4]
  NCORE=args[5]
  message("Input: ", args)
} else if (length(args) == 4) {
  org=args[1]
  startCodon=args[2]
  refTranscriptomeFile=args[3]
  gtfFile=args[4]
  message("Input: ", args)
} else if (length(args) == 0) {
  message("Use default settings")
} else {
  stop("Wrong input? Have you given the correct values for orgamism (scientifi name), start codon, transcript sequences, gene annotation and ncore (optional)?")
}

message("orgamism: ", org)
message("start codon: ", startCodon)
message("Transcriptome file (FASTA): ", refTranscriptomeFile)
message("Gene annotation file (GTF): ", gtfFile)
message("ncore: ", NCORE)

suppressPackageStartupMessages(library(riboSeqR))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))

# Find novel ORFs
# The findCDS function defines (from a fasta file of the transcriptome) all possible coding sequences with a start and stop codon in frame.

message("Loading annotation...")
mouse.txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", organism = org)

#startCodonLlist <- c("ATG", "CTG", "TTG", "GTG")

message("Searching for ORFs with start codon: ", startCodon)
allORFs.transcriptome.orig <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c(startCodon), stopCodon = c("TAG", "TAA", "TGA"))
message("Total ORFs: ", length(allORFs.transcriptome.orig), " (initial search)")
  
# Remove ORFs that is not n*3 (n>1) nucleotides long
allORFs.transcriptome.3n <- allORFs.transcriptome.orig[which(width(allORFs.transcriptome.orig) >= 6 & width(allORFs.transcriptome.orig)%%3 == 0)]
message("Total ORFs: ", length(allORFs.transcriptome.3n), " (after removing regions not n*3 nucleotides long)")
  
allORFs.transcriptome <- allORFs.transcriptome.3n
  
# seqlevels contains uniq seqnames
tx.names <- unlist(lapply(strsplit(seqlevels(allORFs.transcriptome), '[|]'), `[[`, 1))
tx.gene.names <- unlist(lapply(strsplit(as.character((seqnames(allORFs.transcriptome))), '[|]'), `[[`, 6))
tx.biotypes <- unlist(lapply(strsplit(as.character((allORFs.transcriptome)), '[|]'), `[[`, 8))
  
tx.exon <- exonsBy(mouse.txdb, by = "tx", use.names=T)[tx.names]
  
seqlevels(allORFs.transcriptome) <- tx.names
  
message("Running mclapply...")
orfToGenimic <- mclapply(1:length(allORFs.transcriptome), function(i) {
  orf <- allORFs.transcriptome[i]
  tx.name <- runValue(seqnames(orf))
    
  tx.map <- pmapFromTranscripts(orf, tx.exon[names(tx.exon) == tx.name])
    
  #names(tx.map) <- paste(tx.name, tx.gene.names[i], tx.biotypes[i], "ORF", start(orf), end(orf), orf$frame, orf$startCodon, sep="_")
    
  #orf.str <- paste("ORF", start(orf), end(orf), orf$frame, sep="_")
  names(tx.map) <- paste(tx.name, tx.gene.names[i], tx.biotypes[i], start(orf), end(orf), orf$frame, orf$startCodon, sep=":")

  tx.map
}, mc.cores=NCORE)
  
grl <- do.call(c, unlist(orfToGenimic, recursive=FALSE))
unlisted <- unlist(grl)
unlisted.hit <- unlisted[which(mcols(unlisted)$hit == T)]
tx.map.hit.bed <- asBED(split(unlisted.hit, names(unlisted.hit)))
  
message("Exporting to BED...")
orfToBED <- export(tx.map.hit.bed, format = "bed")
write(orfToBED, paste("orfs_",startCodon,".bed", sep=""))
