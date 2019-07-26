### Generate pUTR3 for ORFs
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(riboSeqR))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 6) {
  WDIR=args[1]
  GTF=args[2]
  TAXID=args[3]
  RIBOBAM=args[4]
  refTranscriptomeFile=args[5]
  THREAD=args[6]
  message("Input:")
  message("WDIR:\t", args[1])
  message("GTF:\t", args[2])
  message("TAXID:\t", args[3])
  message("RIBOBAM:\t", args[4])
  message("refTranscriptomeFile:\t", args[5])
  message("THREAD:\t", args[6])
} else if (length(args) == 5) {
  WDIR=args[1]
  GTF=args[2]
  TAXID=args[3]
  RIBOBAM=args[4]
  refTranscriptomeFile=args[5]
  THREAD=8
  message("Input:")
  message("WDIR:\t", args[1])
  message("GTF:\t", args[2])
  message("TAXID:\t", args[3])
  message("RIBOBAM:\t", args[4])
  message("refTranscriptomeFile:\t", args[5])
  } else {
  stop("Input wrong? The args should be [1] working dir [2] GTF [3] TAXID [4] RIBO BAM [5] refTranscriptomeFile [6] Thread")
}

setwd(WDIR)

#refTranscriptomeFile <- "ORFScore_ALL_filtered_expressed_label_regionFiltered_tx.fa"

allORFs.transcriptome.orig.ATG <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
allORFs.transcriptome.orig.CTG <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("CTG"), stopCodon = c("TAG", "TAA", "TGA"))
allORFs.transcriptome.orig.TTG <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("TTG"), stopCodon = c("TAG", "TAA", "TGA"))
allORFs.transcriptome.orig.GTG <- findCDS(fastaFile = refTranscriptomeFile, startCodon = c("GTG"), stopCodon = c("TAG", "TAA", "TGA"))

allORFs.transcriptome <- c(allORFs.transcriptome.orig.ATG, allORFs.transcriptome.orig.CTG, allORFs.transcriptome.orig.TTG, allORFs.transcriptome.orig.GTG)

### Filter ORFs by size
allORFs.transcriptome.filteredBySize <- allORFs.transcriptome[which(width(allORFs.transcriptome) >= 6 & width(allORFs.transcriptome)%%3 == 0)]
allORFs.transcriptome.filteredBysize.sorted <- sort(allORFs.transcriptome.filteredBySize)

### Create ORFId
tx.Ids <- unlist(lapply(strsplit(as.character((allORFs.transcriptome.filteredBysize.sorted)), '[|]'), `[[`, 1))
tx.geneSymbols <- unlist(lapply(strsplit(as.character((allORFs.transcriptome.filteredBysize.sorted)), '[|]'), `[[`, 6))
tx.len <- as.integer(unlist(lapply(strsplit(as.character((allORFs.transcriptome.filteredBysize.sorted)), '[|]'), `[[`, 7)))
tx.biotypes <- unlist(lapply(strsplit(as.character((allORFs.transcriptome.filteredBysize.sorted)), '[|]'), `[[`, 8))
allORFs.transcriptome.filteredBysize.sorted$txId <- tx.Ids
allORFs.transcriptome.filteredBysize.sorted$txLen <- tx.len
allORFs.transcriptome.filteredBysize.sorted$ORFId <- paste(tx.Ids, tx.geneSymbols, tx.biotypes, start(allORFs.transcriptome.filteredBysize.sorted), end(allORFs.transcriptome.filteredBysize.sorted), allORFs.transcriptome.filteredBysize.sorted$frame, allORFs.transcriptome.filteredBysize.sorted$startCodon, sep=":")

#seqlevels(allORFs.transcriptome.filteredBysize.sorted, force=TRUE) <- as.character(unique(seqnames(allORFs.transcriptome.filteredBysize.sorted)))
seqlevels(allORFs.transcriptome.filteredBysize.sorted) <- as.character(unique(seqnames(allORFs.transcriptome.filteredBysize.sorted)))

allORFs.transcriptome.filteredBysize.sorted.split <- split(allORFs.transcriptome.filteredBysize.sorted, seqnames(allORFs.transcriptome.filteredBysize.sorted))
save(allORFs.transcriptome.filteredBysize.sorted.split, file=paste0(refTranscriptomeFile, ".rda"))

#cached_file_names=as.list(dir("../../prediction/", pattern="*.rda$"))
#lapply(cached_file_names,load,.GlobalEnv)
	
#for(i in 1:length(cached_file_names)) {
#    startCodon <- strsplit(strsplit(cached_file_names[i][[1]], "_")[[1]][2], "[.]")[[1]][1]
#    load(paste0("../../prediction/mouse/", cached_file_names[i][[1]]))
#    assign(paste0("allORFs.transcriptome.filteredBysize.sorted.split.", startCodon), allORFs.transcriptome.filteredBysize.sorted.split)
#    rm(allORFs.transcriptome.filteredBysize.sorted.split)
#}

#allORFs.transcriptome.filteredBysize.sorted.split = c(allORFs.transcriptome.filteredBysize.sorted.split.ATG, allORFs.transcriptome.filteredBysize.sorted.split.CTG, allORFs.transcriptome.filteredBysize.sorted.split.GTG, allORFs.transcriptome.filteredBysize.sorted.split.TTG)


#----------------

r <- mclapply(1:length(allORFs.transcriptome.filteredBysize.sorted.split), function(i) {
  GObj <- allORFs.transcriptome.filteredBysize.sorted.split[[i]]
  #message(i)
  #if(length(GObj) > 0) {
  #trLen <- as.integer(strsplit(as.character(runValue(seqnames(GObj[1]))), "[|]")[[1]][7])
  #GObj <- GObj[which(width(GObj) >= 6 & width(GObj)%%3 == 0)]
  
  idx.list <- unlist(lapply(1:length(GObj), function(j) {idx <- which(end(GObj[j]) < start(GObj))[1]}))
  
  idx.na <- which(is.na(idx.list))
  idx.list[is.na(idx.list)] <- 1
  
  pseudo3p.start <- end(GObj) + 1
  pseudo3p.end <- start(GObj[idx.list]) - 1 
  pseudo3p.end[idx.na] <- GObj[1]$txLen
  
  gr <- GRanges(seqnames=seqnames(GObj),
                ranges=IRanges(pseudo3p.start, pseudo3p.end),
                strand="*",
                frame=GObj$frame,
                startCodon=GObj$startCodon,
                stopCodon=GObj$stopCodon,
                context=GObj$context,
                minus3=GObj$minus3,
                plus1=GObj$plus1,
                fromORF=IRanges(ranges(GObj)),
                ORFId=GObj$ORFId,
                txId=GObj$txId
  )
  #}
  #cat("\r",i)
}, mc.cores = THREAD, mc.preschedule=F) 

pUTR3 <- do.call(c, r)

#----------
orfs_filtered <- import("filtered_nonCanonical.bed")
orfs_filtered_txId <- unique(unlist(lapply(strsplit(orfs_filtered$name, '[:]'), `[[`, 1)))

pUTR3.filtered <- pUTR3[pUTR3$ORFId %in% unique(orfs_filtered$name)]
pUTR3.filtered$UTR3Id <- paste(pUTR3.filtered$ORFId,"pUTR3", start(pUTR3.filtered), end(pUTR3.filtered), sep=":")
pUTR3.filtered$UTR3_len <- end(pUTR3.filtered) - start(pUTR3.filtered) + 1

pUTR3.filtered.noUTR.ORFId <- pUTR3.filtered[pUTR3.filtered$UTR3_len <= 0]$ORFId
pUTR3.filtered.withUTR <- pUTR3.filtered[pUTR3.filtered$UTR3_len > 0]
save(pUTR3.filtered.withUTR, file = "pUTR3_info.RData")
#write.table(pUTR3.filtered.withUTR$ORFId, "orf_withUTR_Id.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


txdb <- makeTxDbFromGFF(GTF, format = "gtf", taxonomyId = as.numeric(TAXID))
tx.exon <- exonsBy(txdb, by = "tx", use.names=T)[orfs_filtered_txId]

#seqlevels(allORFs.transcriptome) <- orfs.biotype.filtered.txId

utrToGenomic <- mclapply(1:length(pUTR3.filtered.withUTR), function(i) {
  utr3 <- pUTR3.filtered.withUTR[i]
  tx.map <- pmapFromTranscripts(utr3, tx.exon[which(names(tx.exon) == utr3$txId)])
  names(tx.map) <- utr3$UTR3Id
  
  #orf.str <- paste("ORF", start(orf), end(orf), orf$frame, sep="_")
  #names(tx.map) <- paste(tx.name, tx.gene.names[i], tx.biotypes[i], orf.str, orf$startCodon, sep=":")
  #cat("\r",i)
  tx.map
}, mc.cores = THREAD)

grl <- do.call(c, unlist(utrToGenomic, recursive=FALSE))
unlisted <- unlist(grl)
unlisted.hit <- unlisted[which(mcols(unlisted)$hit == T)]
tx.map.hit.bed <- asBED(split(unlisted.hit, names(unlisted.hit)))
#tx.map.hit.bed$orfName <- unlist(strsplit(tx.map.hit.bed$name, "_pUTR3"))
pUTR3ToBED <- export(tx.map.hit.bed, format = "bed")
write(pUTR3ToBED, file = "pUTR3_BED12.bed")
system("bed12ToBed6 -i pUTR3_BED12.bed > pUTR3_BED6.bed")

write.table(pUTR3.filtered.noUTR.ORFId, "orf_noUTR_Id.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Done!")
