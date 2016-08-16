#!/usr/bin/Rscript

# Author: Fengyuan Hu
# Search ORFs from a given transcript sequences, covert local coordinates to genomic coordinates in BED format
# Input:
# Transcriptome (protein coding + lincRNAs) in Fasta format
# Transcriptome (protein coding + lincRNAs) annotation in GTF format

# Make sure annotation is consistent with sequences, or extract sequences from genome?

suppressPackageStartupMessages(library(riboSeqR))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))

setwd("~/tmp")

# Find novel ORFs
# The findCDS function defines (from a fasta file of the transcriptome) all possible coding sequences with a start and stop codon in frame.
refTranscriptomeFile <- "~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/pc_and_lncRNA_transcripts.fa"

message("Searching for ORFs...")
allORFs.transcriptome.orig <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("ATG", "CTG", "TTG", "GTG"), stopCodon = c("TAG", "TAA", "TGA"))
message("Total ORFs: ", length(allORFs.transcriptome.orig), " (initial search)")

# Remove ORFs that is not n*3 nucleotides long
allORFs.transcriptome.3n <- allORFs.transcriptome.orig[which(width(allORFs.transcriptome.orig)%%3 == 0)]
message("Total ORFs: ", length(allORFs.transcriptome.3n), " (after removing regions not n*3 nucleotides long)")

message("Loading annotation...")
gtfFile <- "~/ref/GENCODE/mouse/M8/annotation/CHR/comprehensive/annotation.gtf"
mouse.txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", organism = "Mus musculus")

###### TESTING ######
## not run ##
#refTranscriptomeFile <- "~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/test/Cxcr4.fa"
#allORFs.transcriptome.orig <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))

#gtfFile <- "~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.chr1.gtf"
#mouse.txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", organism = "Mus musculus")
###### TESTING! ######

allORFs.transcriptome <- allORFs.transcriptome.3n

# seqlevels contains uniq seqnames
tx.names <- unlist(lapply(strsplit(seqlevels(allORFs.transcriptome), '[|]'), `[[`, 1))
tx.exon <- exonsBy(mouse.txdb, by = "tx", use.names=T)[tx.names]

seqlevels(allORFs.transcriptome) <- unlist(lapply(strsplit(seqlevels(allORFs.transcriptome), '[|]'), `[[`, 1))

message("Running mclapply...")
orfToGenimic <- mclapply(1:length(allORFs.transcriptome), function(i) {
  orf <- allORFs.transcriptome[i]
  tx.name <- runValue(seqnames(orf))
  tx.map <- pmapFromTranscripts(orf, tx.exon[names(tx.exon) == tx.name])
  
  names(tx.map) <- paste(tx.name, "ORF", start(orf), end(orf), mcols(orf)$frame, sep="_")
  tx.map
}, mc.cores=4)

grl <- do.call(c, unlist(orfToGenimic, recursive=FALSE))
unlisted <- unlist(grl)
unlisted.hit <- unlisted[which(mcols(unlisted)$hit == T)]
tx.map.hit.bed <- asBED(split(unlisted.hit, names(unlisted.hit)))

message("Exporting to BED...")
orfToBED <- export(tx.map.hit.bed, format = "bed")
write(orfToBED, "orfs.bed")
