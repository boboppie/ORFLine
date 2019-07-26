# To find regions nested in a bigger region

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  WDIR=args[1]
  THREAD=args[2]
  message("Input:")
  message("WDIR:\t", args[1])
  message("THREAD:\t", args[2])
} else if (length(args) == 1) {
  WDIR=args[1]
  THREAD=8
  message("Input: ")
  message("WDIR:\t", args[1])
  } else {
  stop("Input wrong? The args should be [1] working dir [2] Thread")
}

setwd(WDIR)

#orfs_filtered <- import("rrs_filter/ORFScore_ALL_filtered_expressed_label_regionFiltered_rrsWithUTRFiltered_rrsMerged_BED12Plus.bed")
orfs_filtered <- import("BED12Plus.bed")

colnames(mcols(orfs_filtered))[5] <- "ORFScore"
colnames(mcols(orfs_filtered))[16] <- "head_tail_ratio"

#colnames(mcols(orfs.biotype.filtered)) <- c("name","score","itemRgb","ORFScore","p-val","q-val","effective_len","coverage","transcript_id","transcript_biotype","gene_name","gene_biotype","orf_len","start_codon","thick","blocks")
#colnames(mcols(orfs.biotype.filtered)) <- c("name","score","itemRgb","ORFScore","p-val","q-val","effective_len","coverage","transcript_id","transcript_biotype","gene_name","orf_len","start_codon","thick","blocks")

#orfs.biotype.filtered.uniq <- unique(orfs.biotype.filtered)

orfs_filtered$transcript_id <- unlist(lapply(strsplit(orfs_filtered$name, '[:]'), `[[`, 1))
orfs_filtered$orf_end_pos <- unlist(lapply(strsplit(orfs_filtered$name, '[:]'), `[[`, 5))

orfs_filtered$endId <- paste(orfs_filtered$transcript_id, orfs_filtered$orf_end_pos, sep = ":")

write.table(orfs_filtered, file = "withEndId.tsv", row.names = F, col.names = F, quote = F, sep = "\t")

#orfs_filtered.dupEnd <- orfs_filtered[which(duplicated(orfs_filtered$endId))]
#orfs.biotype.filtered.dupEnd[which(orfs.biotype.filtered.dupEnd$chrEnd == "chr10:105574339")]
#split(orfs.biotype.filtered.uniq.dupEnd, orfs.biotype.filtered.uniq.dupEnd$chrEnd)
orfs_filtered_splitByEndId <- split(orfs_filtered, orfs_filtered$endId)

r <- mclapply(1:length(orfs_filtered_splitByEndId), function(i) {
  GObj <- orfs_filtered_splitByEndId[[i]]
  if (length(GObj) == 1) {
    retunrGr <- GObj[1]
  } else {
    # We test start codon first, ATG has the highest priority
    idx.startCodon <- which(endsWith(GObj[]$name, "ATG"))
    if (length(idx.startCodon) == 0) {
      # No ATG, max ORFScore is selected
      retunrGr <- GObj[which.max(GObj[]$ORFScore)]      
    } else {
      if (length(idx.startCodon) == 1) {
        retunrGr <- GObj[idx.startCodon]
      } else {
        # More than 1 ATG, max ORFScore is selected
        retunrGr <- GObj[idx.startCodon][which.max(GObj[idx.startCodon]$ORFScore)]  
      }
    }

    #idx <- which.max(GObj[]$head_tail_ratio)
    #retunrGr <- GObj[idx]
  }
  
  retunrGr
  #cat("\r",i)
}, mc.cores = THREAD) 

nested_filtered <- do.call(c, r)
write.table(nested_filtered, file = "nested_filtered.tsv", row.names = F, col.names = F, quote = F, sep = "\t")

#n_occur <- data.frame(table(orfs.biotype.filtered.uniq$chrEnd))
#n_occur[n_occur$Freq > 1,]$Var1
#orfs.biotype.filtered.uniq[orfs.biotype.filtered.uniq$chrEnd %in% n_occur[n_occur$Freq > 1,]$Var1]

#orfs.biotype.filtered.uniq.noccur <- orfs.biotype.filtered.uniq[orfs.biotype.filtered.uniq$chrEnd %in% n_occur[n_occur$Freq > 1,]$Var1]
#orfs.biotype.filtered.uniq.noccur.split <- split(orfs.biotype.filtered.uniq.noccur, orfs.biotype.filtered.uniq.noccur$chrEnd)

message("Done!")
