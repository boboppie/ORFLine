# visualize ribosome profiling data

library(ggplot2)
library(grid)

# Counts by position of ORFs
setwd("tmp")
df.raw <- read.delim("ENSMUST00000016639.11_ORF_246_833_2-chr19-29360816-29360927.txt", header = F)

pos <- as.factor(c(29360817:29360927))
count <- df.raw$V1
per <- count/sum(count)

df.dummy <- data.frame(pos, count, per)

# Bar plot
p1 <- ggplot(df.dummy, aes(pos, per)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(breaks = c("29360817", "29360927")) +
  theme_bw ()

# Color blocks
per.mat <- matrix(df.dummy$per, nrow = 3)
per.mat.max <- apply(per.mat, 2, max)
per.mat.rep <- t(replicate(3, per.mat.max))
per.rep <- as.vector(per.mat.rep)

per.mat.max.col <- max.col(t(per.mat))
per.mat.max.col.rep <- t(replicate(3, per.mat.max.col))
frame <- as.factor(as.vector(per.mat.max.col.rep))

df.dummy <- data.frame(df.dummy, per.rep, frame)

p2 <- ggplot(df.dummy) + 
  geom_tile(aes(pos, rep("band-per", length(per.rep)), fill = per.rep ,  height = 0.05)) +
  scale_fill_gradient(low="gray90", high="blue") +
  scale_x_discrete(breaks = c("29360817", "29360927")) + theme_minimal() + 
  theme(axis.title.x = element_blank(),axis.text.y = element_text(angle=90))


p3 <- ggplot(df.dummy) + 
  geom_tile(aes(pos, rep("band-frame", length(per.rep)), fill = frame,  height = 0.05)) +
  scale_fill_manual(breaks = c("0", "1", "2", "3"), values = c("#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  scale_x_discrete(breaks = c("29360817", "29360927")) + theme_minimal() + 
  theme(axis.title.x = element_blank(),axis.text.y = element_text(angle=90))

#REF - http://felixfan.github.io/stacking-plots-same-x/

require(gridExtra)
grid.arrange(p3, p2, p1, ncol = 1, heights = c(2,2,5))

ggsave("ORFPlot.png", width = 8, height = 6)

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last"))

