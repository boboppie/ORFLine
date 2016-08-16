# To plot stripchart of sequencing read counts

suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))

setwd("~/Downloads/")

seq.count <- read.table("B_Manuel_sequencing_read_counts_Riboseq_pipeline_v2.tsv", sep="\t", head=T)

seq.count.stack <- melt(seq.count, id = 1:2)
colnames(seq.count.stack)[3:4] <- c("qc_step", "count")

ggplot(seq.count.stack, aes(condition, count, colour=condition)) +
               geom_point(position=position_jitter(width = 0), size=5) + 
               facet_grid(~qc_step) +
               stat_summary(fun.y=mean, aes(ymin=..y.., ymax=..y..), 
               geom='errorbar', width=0.5, color='black', size=1.25) +
               theme_bw()

ggsave("B_Manuel_seq_count_stripchart.pdf", width = 10, height = 6)