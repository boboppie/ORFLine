## ------------------------------------------------------------------------
## TG paper Fig1C
## Input: full coverage, gene-transcripts, problematic regions, target regions, gene list, exon regions

# bash cmd
# $ while read g; do Rscript fig1CandSup.R $g; done <../TG/final/genelist.final.txt
# $ pdfunite *.pdf suppl-63GeneCov.pdf
# $ Rscript fig1CandSup.R NBEAL2 47044566 47048939
## ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
message("Input:")
print(args)

if(length(args) == 1) {
  gene.symbol=args[1]
} else if (length(args) == 3) {
  gene.symbol=args[1]
  window.start=as.integer(args[2])
  window.end=as.integer(args[3])
} else {
  stop("Input wrong?")
}

rm(args)
#gene.symbol="MCFD2"

tx.id <- "ENSMUST00000142893.1"
tx.start <- 128589099
tx.stop <-  128592293 
tx.strand <- "-"
ref <- "mm10"
chr <- "chr1"

if (tx.strand == "-") {complement = T}

itrack <- IdeogramTrack(genome=ref, chromosome=chr, fontcolor="black")
gtrack <- GenomeAxisTrack(fontcolor="black")

ensembl.grcm38.mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                              host="www.ensembl.org",
                              path="/biomart/martservice",
                              dataset="mmusculus_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = ref, chromosome = chr, start = tx.start, end = tx.stop, 
                                    name = "ENSEMBL",  biomart=ensembl.grcm38.mart, transcriptAnnotation = "transcript",
                                    just.group = "above", fontsize=30, fontcolor.group="black", fontsize.title=12)

strack <- SequenceTrack(Mmusculus, chromosome = chr, add53 = TRUE, complement = complement)

track.list = c(itrack, gtrack, biomTrack, strack)
plotTracks(track.list, from = tx.start, to = tx.stop, extend.left=0, extend.right=0,  cex = 1.5)

###########################################################################
## ------------------------------------------------------------------------
setwd(".")
df.header = c("chr", "start", "end", "symbol")

tggenes = read.delim("../TG2-genelist-flank1000.interval", header=F)
names(tggenes) = df.header

tgtargets = read.delim("../TG2_target_mapped_genes.filtered.interval", header=F)
names(tgtargets) = df.header

tgexons = read.delim("../TG2-exon.interval", header=F)
names(tgexons) = df.header

tgpr5p = read.delim("../TG2-problematic.5p.interval", header=F)
names(tgpr5p) = df.header

tgtg = read.delim("../TG2-transcript-gene.info", header=F)
names(tgtg) = c("transcriptid", "symbol")

## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(biomaRt))
# Gviz >= 1.10.11
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))

## ------------------------------------------------------------------------
#genome.hg19 <- BSgenome.Hsapiens.UCSC.hg19
#si=seqinfo(genome.hg19)
#ref=unique(genome(si))
ref="mm10"

## ------------------------------------------------------------------------
geneinfo = tggenes[tggenes$symbol==gene.symbol,]
chr = paste("chr", as.character(unique(geneinfo$chr)), sep="")

# window.start and window.end are the min(target start) and max(target end)
# GFI1B doesn't have target regions???
genetargets = tgtargets[tgtargets$symbol==gene.symbol,]
if (nrow(genetargets) > 0) {
  if (!exists(as.character(substitute(window.start))) && !exists(as.character(substitute(window.end)))) {
    window.start=min(genetargets$start)
    window.end=max(genetargets$end)
  }
  
  gr.targets <- GRanges(seqnames=Rle(c(chr), c(nrow(genetargets))),
                        ranges=IRanges(genetargets$start, genetargets$end),
                        strand="*")
} else {
  if (!exists(as.character(substitute(window.start))) && !exists(as.character(substitute(window.end)))) {
    window.start=geneinfo$start
    window.end=geneinfo$end
  }
}

genepr5p = tgpr5p[tgpr5p$symbol==gene.symbol,]
if(nrow(genepr5p) > 0) {
  gr.pr5p <- GRanges(seqnames=Rle(c(chr), c(nrow(genepr5p))),
                     ranges=IRanges(genepr5p$start, genepr5p$end),
                     strand="*")
}

## ------------------------------------------------------------------------

gene.cov.path=paste("../merge/gene/", gene.symbol, ".cov", sep="")
gene.no_col <- max(count.fields(gene.cov.path, sep = " "))
gene.raw <- read.table(gene.cov.path, sep=" ", fill=TRUE, header=T)
gene.raw[is.na(gene.raw)] <- 0

gene.raw.noloc <- data.frame(gene.raw[2:ncol(gene.raw)])
gene.raw.loc <- data.frame(gene.raw[1])
#colnames(exon.raw.loc) <- c("loc")
gene.raw.noloc.m <- data.matrix(gene.raw.noloc)
gene.cov.median <- rowMedians(gene.raw.noloc.m)
gene.cov.lower <- apply(gene.raw.noloc, 1, quantile, probs = c(0.05),  na.rm = TRUE)
gene.cov.upper <- apply(gene.raw.noloc, 1, quantile, probs = c(0.95),  na.rm = TRUE)

gene.cov.summary <- data.frame(gene.raw.loc, gene.cov.lower, gene.cov.median, gene.cov.upper)
colnames(gene.cov.summary) <- c("loc", "lower", "median", "upper")

# gr means GRanges
df.gr.gene.cov.raw = cbind(data.frame(chr=chr, start=gene.raw.loc$loc, end=gene.raw.loc$loc,
                                      strand="*"), gene.raw.noloc)
gr.gene.cov.raw <- makeGRangesFromDataFrame(df.gr.gene.cov.raw, keep.extra.columns=TRUE)

df.gr.gene.cov.summary = cbind(data.frame(chr=chr, start=gene.cov.summary$loc, end=gene.cov.summary$loc,
                                          strand="*"), gene.cov.summary[,2:4])
gr.gene.cov.summary <- as(df.gr.gene.cov.summary, "GRanges")
dTrack.gene.summary <- DataTrack(gr.gene.cov.summary,
                                 name = "Coverage",
                                 type = c("l", "g"),
                                 groups = c("Lower (5th percentile)", "Median (50th percentile)", "Upper (95th percentile)"),
                                 col=c("red", "#0099FF", "darkgreen"),
                                 #cex = 0.5, 
                                 legend = TRUE,
                                 #box.legend=T,
                                 #pch=20,
                                 cex.legend=0.8,
                                 fontsize.legend=30,
                                 #fontsize.title=12,
                                 
                                 fontcolor.legend="black",
                                 cex.axis=1.7
)

## ------------------------------------------------------------------------
ensembl.grch37.mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                              host="grch37.ensembl.org",
                              path="/biomart/martservice",
                              dataset="hsapiens_gene_ensembl")
fm <- Gviz:::.getBMFeatureMap()
fm["symbol"] <- "external_gene_name"

## ------------------------------------------------------------------------

track.list = list(dTrack.gene.summary)

if (exists(as.character(substitute(gr.pr5p)))) {
  atrack.pr5p <- AnnotationTrack(gr.pr5p, name = "Prob. regions", fill="red", stacking = "dense", size=0.65)
  track.list = c(track.list, atrack.pr5p)
}

if(exists(as.character(substitute(gr.targets)))) {
  atrack.target <- AnnotationTrack(gr.targets, name = "Targets",  stacking = "dense")
  track.list = c(track.list, atrack.target)
}

biomTrack <- BiomartGeneRegionTrack(genome=ref, chromosome=chr, start=window.start, end=window.end, name=gene.symbol,
                                    featureMap=fm, filters=list(ensembl_transcript_id=tgtg[tgtg$symbol==gene.symbol,]$transcriptid),
                                    biomart=ensembl.grch37.mart, transcriptAnnotation = "transcript",
                                    just.group = "above", fontsize=30, fontcolor.group="black", fontsize.title=12)
track.list = c(track.list, biomTrack)

itrack <- IdeogramTrack(genome=ref, chromosome=chr, fontsize=30, fontcolor="black")
gtrack <- GenomeAxisTrack(fontsize=20, fontcolor="black")

track.list = c(track.list, gtrack, itrack)

#plotTracks(list(itrack, gtrack, dTrack.gene.summary),  background.panel = "#FFFEDB", background.title = "darkblue", cex.title=1.5)

#biomTrack <- BiomartGeneRegionTrack(genome=ref, chromosome=chr, start=window.start, end=window.end, name="ENSEMBL GRCh37 Gene Models",
#                                    featureMap=fm, biomart=ensembl.grch37.mart, transcriptAnnotation = "symbol")

## ------------------------------------------------------------------------
#w = 26
h = 14

# plot width is propotional to gene length, set a factor, p.fac=1500 based on length of ANXA5 (~=30000) / width (20)
p.fac=1500
w=(window.end-window.start+2000)/p.fac
if (w<26) {w=26}

pdf(paste(gene.symbol, ".pdf", sep=""), width=w, height=h)
plotTracks(track.list, frame=F, from = window.start-1000, to = window.end+1000,
           #background.panel = "#FFFEDB", background.title = "darkblue", 
           #cex=5,
           col.axis="black", col.title="black", background.title="lightgray",
           cex.title=1.5)
dev.off()