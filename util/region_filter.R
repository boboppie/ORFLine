#!/usr/bin/Rscript

### Author: Fengyuan Hu
### R code to filter out duplcated, overlapped regions

suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library(GenomicAlignments))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 4) {
  WDIR=args[1]
  CDSGTF=args[2]
  BAM=args[3]
  THREAD=args[4]
  message("Input:")
  message("WDIR:\t", args[1])
  message("CDSGTF:\t", args[2])
  message("BAM:\t", args[3])
  message("THREAD:\t", args[4])
} else if (length(args) == 3) {
  WDIR=args[1]
  CDSGTF=args[2]
  BAM=args[3]
  THREAD=8
  message("Input:")
  message("WDIR:\t", args[1])
  message("CDSGTF:\t", args[2])
  message("BAM:\t", args[3])
  message("THREAD:\t", 8)
  } else {
  stop("Input wrong? The args should be [1] working dir [2] CDS gtf [3] BAM [4] Thread")
}

setwd(WDIR)

cds <- import(CDSGTF)
cds.uniq <- unique(cds)
message("GTF loaded.")

ribo_bam <- readGAlignments(BAM, use.names = TRUE)
message("BAM loaded.")

cds.uniq.olap <- findOverlaps(cds.uniq, ribo_bam)
#cds.uniq$ribo_hits <- tabulate(queryHits(cds.uniq.olap), length(cds.uniq))
#cds.uniq$ribo_hits <- countQueryHits(cds.uniq.olap)
#cds.uniq.hasHits <- cds.uniq[cds.uniq$ribo_hits>0]

# Using BED12 regions will take reads aligned to intronic regions into account, inflated the number
# use bed12ToBed6 -i biotype_filtered_ATG_BED12.bed >biotype_filtered_ATG_BED6.bed

#orfs.olap <- findOverlaps(orfs, ribo_bam)
#orfs$ribo_hits <- tabulate(queryHits(orfs.olap), length(orfs))

#codon <- "ATG"
orfs  <- import("ORFScore_ALL_filtered_expressed_label_BED12Plus.bed")
orfs_bed6 <- import("ORFScore_ALL_filtered_expressed_label_BED6.bed")
message("BED loaded.")

orfs_bed6_splited <-split(orfs_bed6, orfs_bed6$name)
orfs_bed6_splited.olap <- findOverlaps(orfs_bed6_splited, ribo_bam)
orfs_bed6_splited.df <- data.frame(names(orfs_bed6_splited))
colnames(orfs_bed6_splited.df) <- "name"
# This gives unique hits
orfs_bed6_splited.df$ribo_hits_nonDup <- countQueryHits(orfs_bed6_splited.olap)
orfs_bed6_splited.df[] <- lapply(orfs_bed6_splited.df, as.character) # factor to string

# orfs_bed6 has to have the same order as orfs_bed12
orfs$ribo_hits_nonDup <- as.integer(orfs_bed6_splited.df$ribo_hits_nonDup[order(match(orfs_bed6_splited.df$name,orfs$name))]) 
message("ORF overlap finished.")

#orfs_bed6_olap <- findOverlaps(orfs_bed6, ribo_bam)

orfs_cds_olap <- findOverlaps(orfs_bed6, cds.uniq)

cds.uniq.olap.gr <- cds.uniq[subjectHits(orfs_cds_olap)]
elementMetadata(cds.uniq.olap.gr)$overlapped_ORFId <- orfs_bed6[queryHits(orfs_cds_olap)]$name
#cds.uniq.olap.gr.uniq <- unique(cds.uniq.olap.gr) # do not do unique, ORFs overlapping with the same exon will be removed
cds.uniq.olap.gr_splited <- split(cds.uniq.olap.gr, cds.uniq.olap.gr$overlapped_ORFId)
cds.uniq.olap.gr_splited.olap <- findOverlaps(cds.uniq.olap.gr_splited, ribo_bam)
cds.uniq.olap.gr_splited.df <- data.frame(names(cds.uniq.olap.gr_splited))
colnames(cds.uniq.olap.gr_splited.df) <- "name"

# This gives unique hits
cds.uniq.olap.gr_splited.df$ribo_hits_nonDup <- countQueryHits(cds.uniq.olap.gr_splited.olap)
cds.uniq.olap.gr_splited.df[] <- lapply(cds.uniq.olap.gr_splited.df, as.character) # factor to string
cds.uniq.olap.gr_splited.df$ribo_hits_nonDup <- as.numeric(cds.uniq.olap.gr_splited.df$ribo_hits_nonDup)
message("CDS overlap finished.")

message("To run in parallel...")

# Calculate ORF/CDS overlap ratio
# 

orf_cds_ratio <- mclapply(1:length(orfs), function(i) {
  #sum(cds.uniq[subjectHits(orfs_cds_olap)[which(queryHits(orfs_cds_olap) == i)]]$ribo_hits)/orfs[i]$ribo_hits_nonDup
  # This give us CDS idx
  #p <- subjectHits(orfs_cds_olap)[queryHits(orfs_cds_olap) %in% which(orfs_bed6$name==orfs[i]$name)]
  cds_hits <- cds.uniq.olap.gr_splited.df[cds.uniq.olap.gr_splited.df$name == orfs[i]$name,]$ribo_hits_nonDup  
  #p <- p[!is.na(p)]
  if (length(cds_hits) == 0) {
    r <- NA
    #message(orfs[i]$name)
  } else {
    #r <- sum(cds.uniq[p]$ribo_hits)/orfs[i]$ribo_hits_nonDup
    # This gives unique counts for CDS 
    #r <- length(unique(subjectHits(cds.uniq.olap[queryHits(cds.uniq.olap) %in% p])))/orfs[i]$ribo_hits_nonDup # slower than the following code
    #r <- length(unique(subjectHits(cds.uniq.olap)[queryHits(cds.uniq.olap) %in% p]))/orfs[i]$ribo_hits_nonDup
    r <- cds_hits/orfs[i]$ribo_hits_nonDup
  }
}, mc.cores=THREAD)
message("Parallel run finished")

orfs$orf_cds_ratio <- unlist(orf_cds_ratio)

write.table(orfs, file = "orfs_readRatio.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

colnames(mcols(orfs))[4] <- "orf_label"

orfs_canonical <- orfs[which(orfs$orf_label == "canonical" & orfs$orf_cds_ratio == 1)]
write.table(orfs_canonical, file = "orfs_readRatio_canonical.tsv", sep = "\t", row.names = F, col.names = F, quote = F)

#orfs[which(orfs$NA..12 %in% c("canonical_truncated", "within", "canonical", "canonical_extended") & orfs$orf_cds_ratio==1)]

#orfs_filtered <- orfs[which(orfs$orf_cds_ratio <0.8 | is.na(orfs$orf_cds_ratio))]
#export(orfs_filtered, "ALL_filtered.bed", "bed")

# To improve, check if orf overlaps with 5' or 3' UTR (possible true cases), otherwise intronic regions (FPs)
#orfs_non_pc <- orfs[orfs$transcript_biotype != "protein_coding"]
#orfs_non_pc[which(orfs_non_pc$orf_cds_ratio >0.8)]

message("All done!")