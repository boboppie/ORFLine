#!/usr/bin/python

# DESCRIPTION: Takes as input a BAM alignment file of ribosome profiling reads,
# gene or transcript abundance estimates from Cufflinks (mRNA-Seq or ribo-Seq)
# and a GTF annotation file (UCSC genes.gtf or Ensembl have been tested. From
# these inputs, the normalized read coverage over each transcript is computed
# and the spectral coherence relative to an enriched signal over a single
# reading frame is calculated (eg. 1:0:0, repeating). Additional implementations
# to calculate the FLOSS read distribution and FLOSS score (Ingolia, 2014), and
# ORFScore (Bazzini, 2014) are included. This version is able to handle one
# sample (BAM) at a time, however an additional script (intersect_experiments.py
# is included to merge multiple experiments into a single analysis.

# If you use SPECtre as part of your analysis, please cite the following:
#

# This script was written and tested using Python v2.7.8 on an OS X v10.9.4
# workstation and on a server running RHEL r6.4.

# DEPENDENCIES:
# bx-python:	https://bitbucket.org/james_taylor/bx-python/wiki/Home
# pyfasta:		https://pypi.python.org/pypi/pyfasta/
# rpy:			http://rpy.sourceforge.net/
# ROCR (in R)	https://rocr.bioinf.mpi-sb.mpg.de/

from plastid import Transcript, BED_Reader, BAMGenomeArray, FivePrimeMapFactory, VariableFivePrimeMapFactory, SizeFilterFactory
import os
import numpy

os.chdir('/Users/huf/')

orfBedFile = "tmp/plastid/data/orfs_Cxcr4_test.bed"
riboAlignmentFile_B_resting_chr1_29nt = "tmp/plastid/data/B_resting_chr1_29nt_12offset.bam"
riboAlignmentFile_B_resting_chr1 = "tmp/plastid/data/B_resting_chr1.bam"
psiteFile = "out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-18_10-44-28/plastid/psite/merged_q255_star_genome_p_offsets.txt"

orfs = list(BED_Reader(orfBedFile ,return_type=Transcript))
alignments = BAMGenomeArray(riboAlignmentFile_B_resting_chr1)
#alignments.set_mapping(FivePrimeMapFactory(offset=12))

maprule = VariableFivePrimeMapFactory.from_file(open(psiteFile))
alignments.set_mapping(maprule)

from plastid import SizeFilterFactory
size_filter = SizeFilterFactory(min=29,max=35)
alignments.add_filter("size",size_filter)

# create a holder for phasing
phasing = numpy.zeros(3)
#phasing_orf = numpy.zeros(3)

# start codons are hyper-phased; stop codons can have differnt
# phasing or even be de-phased depending on experimental protocol
# so, we'll ignore 5 codons after the start, and 5 before the stop
#codon_buffer = 5*3
codon_buffer = 0

for my_orf in orfs:
	cds = my_orf.get_cds()
	# if transcript is coding
	 if len(cds) > 0:
	 	try:
	 		# get numpy.ndarray of counts in coding region
	 		# for reverse-strand features, counts are reversed relative to genomic coordinates
	 		#counts = cds.get_counts(alignments)[codon_buffer:-codon_buffer]
	 		counts = cds.get_counts(alignments)

	 		# reshape to Nx3, where N = number of codons
	 		counts = counts.reshape((len(counts)/3,3))

	 		# sum over codon positions to get a 3-vector,
	 		# and add to data holder
	 		
	 		#phasing_orf += counts.sum(0).astype(float) / counts.sum(0).sum()

	 		phasing += counts.sum(0)

	 	except: # raise exception if coding region is not n*3 nucleotides long
	 		print("Length (%s nt) of CDS for `%s` contains partial codons. Frameshift?" % (len(counts),my_transcript.get_name()))

# compute fraction of phased reads
#phasing_proportions = phasing.astype(float) / phasing.sum()
#phasing_proportions