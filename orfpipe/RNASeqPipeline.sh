#!/bin/bash

# Parameters
CELLTYPE=B
CONDITION=Resting
AUTHOR=manuel

SAMPLEID=WT1
RAWDATAFILENAME=${CELLTYPE}_rna-seq_${AUTHOR}_${CONDITION}_${SAMPLEID}

THREADS=4

REFGENOME=~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome
SPLICESITESFILE=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.splicesites.txt
PROTEINCODINGGTF=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf
PHASINGGTF=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/plastid.phasing.transcript.gtf


RAWDATAPATH=~/data/orf-discovery/${CELLTYPE}/rna-seq/${AUTHOR}/${CONDITION}/fastq
OUTPATH=~/out/orf-discovery/${CELLTYPE}/rna-seq/${AUTHOR}/${CONDITION}

TIMESTAMP=$(date +"%Y-%b-%d_%H-%M-%S")
OUTPATH=$OUTPATH/$TIMESTAMP

mkdir -p $OUTPATH/cutadapt
mkdir -p $OUTPATH/fastqc
mkdir -p $OUTPATH/hisat2
mkdir -p $OUTPATH/fastqc/hisat2/transcriptome

fastqc -t $THREADS -o $OUTPATH/fastqc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz

# Trim Adapter
# Base "A" is added before normal adapter sequence (GATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
# REF http://cutadapt.readthedocs.io/en/stable/guide.html#warning-about-incomplete-adapter-sequences
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -f fastq -e 0.05 -O 3 --quality-base=33 --trim-n -m 60 --max-n=0.05 -o $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq <(gzip -dc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz) > $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.log

gzip $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq

fastqc -t $THREADS -o $OUTPATH/fastqc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz

# HISAT2 Genome Alignment
hisat2 -q -p $THREADS -x $REFGENOME \
       --known-splicesite-infile $SPLICESITESFILE \
       --seed 23 -U $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz \
       --dta --dta-cufflinks \
       -S $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.log

samtools view -bS $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools index $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools idxstats $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.totalreads

samtools view -q 60 -b $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam
samtools index $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam
samtools idxstats $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam.idxstats
cut -f3 $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam.totalreads

fastqc -t $THREADS -o $OUTPATH/fastqc/hisat2/transcriptome -f bam $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam

rm $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam

# StringTie transcriptome assembly
mkdir -p $OUTPATH/stringtie/${RAWDATAFILENAME}

stringtie $OUTPATH//hisat2/${RAWDATAFILENAME}_transcriptome_q60.sorted.bam \
          -p $THREADS -m 30 -G $PROTEINCODINGGTF \
          -o $OUTPATH/stringtie/${RAWDATAFILENAME}/assembly.gtf \
          -A $OUTPATH/stringtie/${RAWDATAFILENAME}/gene_abund.tab \
          -C $OUTPATH/stringtie/${RAWDATAFILENAME}/cov_refs.gtf

stringtie --merge -G $PROTEINCODINGGTF -o $OUTPATH/stringtie/merged.gtf `find $OUTPATH/stringtie -name assembly.gtf`
awk '{if($3=="transcript")print}' merged.gtf >merged.transcript.gtf

# Extracting transcript sequences using gffread or Galaxy (run gtft2ff3.pl)
# http://cole-trapnell-lab.github.io/cufflinks/file_formats/
# StringTie will output strand information, but what does dot mean? Strand can't be determined?
