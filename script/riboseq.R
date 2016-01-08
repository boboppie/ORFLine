suppressMessages(library(riboSeqR))

setwd("~/tmp/ribo-seq_test/data/manuel/se/bowtie")

ribofile <- "transcriptome_mapped.bowtie"
riboDat <- readRibodata(ribofile, 
                        replicates = "wt1", 
                        columns = c(strand = 1, seqname = 2, start = 3, sequence = 4))
lengthDist(riboDat)

# Use Ensembl transcript_biotype:protein_coding
transcriptomeFasta <- "../../../../ref/ensembl/mouse/transcriptome/Mus_musculus.GRCm38.cdna.protein_coding_transcripts.fa"

# How does this step work??? What about reverse strand?
# fastaCDS <- findCDS(fastaFile = transcriptomeFasta, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
# CDS information from transcript_coding_union.R
CDS <- transcripts.merged.gr
CDS$frame <- (start(CDS)-1)%%3

len <- c(20:40)
fCs <- frameCounting(riboDat, CDS, lengths = len)
fS <- readingFrame(rC = fCs, lengths = len)
fS
plotFS(fS, lengths = len)

