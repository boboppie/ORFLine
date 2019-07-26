### Generate pUTR3 for ORFs
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicAlignments))
#suppressPackageStartupMessages(library(riboSeqR))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  WDIR=args[1]
  RNABAMDIR=args[2]
  message("Input:")
  message("WDIR:\t", args[1])
  message("RNABAMDIR:\t", args[2])
} else if (length(args) == 3) {
  WDIR=args[1]
  RNABAMDIR=args[2]
  CELLLINE=args[3]
  message("Input:")
  message("WDIR:\t", args[1])
  message("RNABAMDIR:\t", args[2])
  message("CELLLINE:\t", args[3])
} else {
  stop("Input wrong? The args should be [1] working dir [2] RNA BAM dir")
}

setwd(WDIR)

#------------------- RRS implementation ---------------------
# input: merged Ribo-seq BAM, merged RNA-seq BAM, ORF, UTR3

message("Loading RNA-seq BAMs ...")

#if (exists("CELLLINE")) {
#    rna_bam_files <- as.list(dir(RNABAMDIR, pattern=paste0('^',CELLLINE,'.*bam$')))
#} else {
#    rna_bam_files <- as.list(dir(RNABAMDIR, pattern="merged_q255.bam$"))
#}

#assign("rna_bam", readGAlignments(paste0(RNABAMDIR, rna_bam_files[1][[1]]), use.names = TRUE))

if (exists("CELLLINE")) {
    rna_bam_files <- as.list(dir(RNABAMDIR, pattern=paste0('^',CELLLINE,'.*bam$')))
} else {
    rna_bam_files <- as.list(dir(RNABAMDIR, pattern="*_q255.sorted.bam$"))
}

# Loading all RNA-seq BAMs for the experiment (2-5 depends on the experiments)
for (i in 1:length(rna_bam_files)) {
    assign(paste0("rna_bam_", i),  readGAlignments(paste0(RNABAMDIR, rna_bam_files[i][[1]]), use.names = TRUE))
}

message("RNA-seq BAM loaded...")

load("pUTR3_info.RData")
ribo_orf_sum_merged <- read.delim("Ribo_ALL_ORF_sum.tsv", header=F)
ribo_utr3_sum_merged <- read.delim("Ribo_ALL_UTR3_sum.tsv", header=F)

orfs <- import("filtered_BED6.bed")
ORFId_canonical <- read.table("ORFId_canonical.txt", header=F)
orfs.withUTR <- orfs[orfs$name %in% c(pUTR3.filtered.withUTR$ORFId, as.character(ORFId_canonical$V1))]
orfs.withUTR.split <-split(orfs.withUTR, orfs.withUTR$name)
orfs.withUTR.split.df <- data.frame(names(orfs.withUTR.split))
colnames(orfs.withUTR.split.df) <- "ORFId"
orfs.withUTR.split.df[] <- lapply(orfs.withUTR.split.df, as.character)

utr3_nonCanonical <- import("pUTR3_BED6.bed")
utr3_canonical <- import("UTR3_canonical_BED6.bed")
utr3 <- c(utr3_nonCanonical, utr3_canonical)
utr3.split <-split(utr3, utr3$name)
utr3.split.df <- data.frame(names(utr3.split))
colnames(utr3.split.df) <- "UTR3Id"
utr3.split.df[] <- lapply(utr3.split.df, as.character)

# The right way to implement RRS, computed the 3'UTR scores as any read overlapping the 3'UTR 
# (excluding reads that also overlap the coding region which will artificially inflate its estimate)
strand.value <- as.character(as.data.frame(strand(utr3))$value)
strand.forw.idx <- which(strand.value == "+")
strand.rev.idx <- which(strand.value == "-")
utr3.bound <- utr3
end(utr3.bound[strand.forw.idx]) <- start(utr3.bound[strand.forw.idx])
start(utr3.bound[strand.rev.idx]) <- end(utr3.bound[strand.rev.idx])
utr3.bound.split <-split(utr3.bound, utr3.bound$name)

orfs.withUTR.split.df$orf_hits_nonDup_ribo_merged <- ribo_orf_sum_merged$V2[trimws(as.character(ribo_orf_sum_merged$V1)) %in% orfs.withUTR$name]

# Translation efficiency need the merged RPF count
df.TE <- orfs.withUTR.split.df

# Add ORFId to utr3.split.df dataframe utr3.split.df
utr3.split.df$ORFId <- unlist(lapply(strsplit(as.character(utr3.split.df$UTR3Id), ":pUTR3"), `[[`, 1))
# Add utr3_hits_nonDup_ribo_merged to datafame utr3.split.df
utr3.split.df$utr3_hits_nonDup_ribo_merged <- ribo_utr3_sum_merged$V2[trimws(as.character(ribo_utr3_sum_merged$V1)) %in% utr3$name]

# Loop through RNA bams to count reads around UTR3
for (i in 1:length(rna_bam_files)) {
    orfs.withUTR.split.df[paste0("orf_hits_nonDup_rna_", i)] <- countQueryHits(findOverlaps(orfs.withUTR.split, get(paste0("rna_bam_", i))))
    utr3.split.df[paste0("utr3_hits_nonDup_rna_", i)] <- countQueryHits(findOverlaps(utr3.split, get(paste0("rna_bam_", i))))
    #utr3.split.df[paste0("utr3_hits_nonDup_rna_all_", i)] <- countQueryHits(findOverlaps(utr3.split, get(paste0("rna_bam_", i))))
    #utr3.split.df[paste0("utr3_hits_nonDup_rna_bound_", i)] <- countQueryHits(findOverlaps(utr3.bound.split, get(paste0("rna_bam_", i))))
    #utr3.split.df[paste0("utr3_hits_nonDup_rna_", i)] <- utr3.split.df[paste0("utr3_hits_nonDup_rna_all_", i)] - utr3.split.df[paste0("utr3_hits_nonDup_rna_bound_", i)]
}

# Merged ORF RNA hits of all samples
orfs.withUTR.split.df$orf_hits_nonDup_rna_merged <- rowSums(orfs.withUTR.split.df[, grep('^orf_hits_nonDup_rna_\\d+$', colnames(orfs.withUTR.split.df)), drop=F])

# Merged UTR3 RNA hits of all samples
utr3.split.df$utr3_hits_nonDup_rna_merged <- rowSums(utr3.split.df[, grep('^utr3_hits_nonDup_rna_\\d+$', colnames(utr3.split.df)), drop=F])

# Merge two dataframes based on ORFId, slim down utr3.split.df to 4 cols, orfs.withUTR.split.df to 3 cols
df.merged <- merge(orfs.withUTR.split.df[, c(1:2, ncol(orfs.withUTR.split.df))], utr3.split.df[,c(1:3, ncol(utr3.split.df))], by="ORFId")

df.merged$orf_len <- unlist(lapply(strsplit(df.merged$ORFId,":"), function(x) {as.numeric(x[5])-as.numeric(x[4])+1}))
df.merged$utr3_len <- read.delim("UTR3_length.tsv", header=F)$V2
df.merged$utr3_orf_len_ratio <- df.merged$utr3_len/df.merged$orf_len

# Calculate RRS
ribo_utr3 <- replace(df.merged$utr3_hits_nonDup_ribo_merged, df.merged$utr3_hits_nonDup_ribo_merged==0, 0.001)
df.merged$ribo_ratio_merged <- (df.merged$orf_hits_nonDup_ribo_merged/ribo_utr3) * df.merged$utr3_orf_len_ratio

rna_utr3 <- replace(df.merged$utr3_hits_nonDup_rna_merged, df.merged$utr3_hits_nonDup_rna_merged==0, 0.001)
df.merged$rna_ratio_merged <- df.merged$orf_hits_nonDup_rna_merged/rna_utr3 * df.merged$utr3_orf_len_ratio

rna_ratio <- replace(df.merged$rna_ratio_merged, df.merged$rna_ratio_merged==0, 0.001)
df.merged$RRS_merged <- df.merged$ribo_ratio_merged/rna_ratio

write.table(df.merged, "RRS_out_full.tsv", row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
write.table(df.merged[,c("ORFId", "ribo_ratio_merged", "rna_ratio_merged", "RRS_merged")], "RRS_out.tsv", row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")

#==============================================================
#
# Translation efficiency is defined as:
# 
# (density of RPF within ORF) / (RNA expression of ORFs transcript)
# 
# Density of RPF is the RPKM value
#
# TE = RP_rpkm / RNA-seq_rpkm (raw, or take log2)
# 
# Note: RPKM (SE) and FPKM (PE) are the same thing 
#
# Calculate RPKM:
#
# RPKM =  numReads/(geneLength/1000 * totalNumReads/1,000,000)
#
# numReads - number of reads mapped to a gene sequence
# geneLength - length of the gene sequence
# totalNumReads - total number of mapped reads of a sample
#
#===============================================================

# smORF length
orf.len <- read.delim("filtered_smORFs_length.tsv", header=F)
colnames(orf.len) <- c("ORFId", "orf_len")
orf.len.sub <- orf.len[orf.len$ORFId %in% orfs.withUTR.split.df$ORFId,]
orf.len.sub.sorted <- orf.len.sub[with(orf.len.sub, order(ORFId)), ]

df.TE$orf_len <- orf.len.sub.sorted$orf_len

rp.totalMappedReads <- read.delim("rp_totalMappedReads_merged.tsv", header=F)
colnames(rp.totalMappedReads) <- "rp_total_merged"

df.TE$rp_rpkm_merged <- df.TE$orf_hits_nonDup_ribo_merged/((df.TE$orf_len/1000) * ((rp.totalMappedReads$rp_total_merged/1000000)[[1]]))

rna_fpkm_merged <- read.delim("rna_fpkm_merged.tsv", header=F)
colnames(rna_fpkm_merged) <- c("TxId", "rna_fpkm_merged")

df.TE$TxId <- unlist(lapply(strsplit(as.character(df.TE$ORFId), ":"), `[[`, 1))

df.TE <- merge(df.TE, rna_fpkm_merged, by="TxId")

df.TE$te_merged <- df.TE$rp_rpkm_merged/df.TE$rna_fpkm_merged
df.TE$te_merged_log2 <- log2(df.TE$te_merged)

write.table(df.TE, "TE_smORF_full.tsv", row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
write.table(df.TE[,c("TxId", "ORFId", "te_merged", "te_merged_log2")], "TE_smORF_teMerged.tsv", row.names = FALSE, col.names = F, quote = FALSE, sep = "\t")


#============================================================================================
#
# CDS TE of uORF
#
# uORFs regulate downstream CDS, different from smORF
# 
# the calculation is - 
# TE = RP_CDS_rpkm / RNA-seq_exon_rpkm
#
#============================================================================================

# import RPF counts for CDS
rp_CDS_PC_merged <- read.delim("CDS_PC_RPF_merged.tsv", header=F)
colnames(rp_CDS_PC_merged) <- c("TxId", "orf_hits_nonDup_ribo_merged")

# import CDS length
df.TE_CDS <- read.delim("CDS_PC_length.tsv", header=F)
colnames(df.TE_CDS) <- c("TxId", "cds_len")

df.TE_CDS$orf_hits_nonDup_ribo_merged <- rp_CDS_PC_merged$orf_hits_nonDup_ribo_merged

df.TE_CDS$rp_rpkm_merged <- df.TE_CDS$orf_hits_nonDup_ribo_merged/((df.TE_CDS$cds_len/1000) * ((rp.totalMappedReads$rp_total_merged/1000000)[[1]]))

df.TE_CDS <- merge(df.TE_CDS, rna_fpkm_merged, by="TxId")
df.TE_CDS$te_merged <- df.TE_CDS$rp_rpkm_merged/df.TE_CDS$rna_fpkm_merged
df.TE_CDS$te_merged_log2 <- log2(df.TE_CDS$te_merged)


write.table(df.TE_CDS, "TE_CDS_PC.tsv", row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
write.table(df.TE_CDS[,c("TxId", "te_merged", "te_merged_log2")], "TE_CDS_PC_teMerged.tsv", row.names = FALSE, col.names = F, quote = FALSE, sep = "\t")

message("Done!")
