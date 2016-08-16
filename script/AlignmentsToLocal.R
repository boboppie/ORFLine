#!/usr/bin/Rscript

# Author: Fengyuan Hu
# Map genomic alignments to transcipt local coordinates

suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))

# Example:
alignments <- GAlignments("chr1", 10L, "11M", strand("*"), names="read_A")
x <- GRanges("chr1", IRanges(c(12, 12), width=c(6, 20)))
mapToAlignments(x, alignments)