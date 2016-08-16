#!/bin/bash

RAWDATAPATH=~/data/orf-discovery/B/ribo-seq/manuel/Resting/fastq
OUTPATH=~/out/orf-discovery/B/ribo-seq/manuel/Resting

#RAWDATAFILENAME=B_ribo-seq_manuel_LPS_WT1
#RAWDATAPATH=~/data/orf-discovery/B/ribo-seq/manuel/LPS/fastq
#OUTPATH=~/out/orf-discovery/B/ribo-seq/manuel/LPS

#REFRIBOFILTER=~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/bt2/ribo-seq_reads-filter-extra-plus-human
REFRIBOFILTER=~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/bt2/ribo-seq_reads-filter
REFGENOME=/Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome
REFTRANSCRIPTOME=/Users/huf/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/pc_and_lncRNA_transcripts
SPLICESITESFILE=/Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.splicesites.txt
PROTEINCODINGGTF=/Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf
PHASINGGTF=/Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/plastid.phasing.transcript.gtf
PHASINGBED=/Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/plastid.phasing.cds.gtf.bed

TIMESTAMP=$(date +"%Y-%b-%d_%H-%M-%S")
OUTPATH=$OUTPATH/$TIMESTAMP
mkdir -p $OUTPATH

mkdir -p $OUTPATH/cutadapt/log
mkdir -p $OUTPATH/fastqc
mkdir -p $OUTPATH/bowtie2/log
mkdir -p $OUTPATH/hisat2/log
mkdir -p $OUTPATH/fastqc/hisat2/transcriptome
mkdir -p $OUTPATH/fastqc/bowtie2/transcriptome
mkdir -p $OUTPATH/plastid/psite
mkdir -p $OUTPATH/plastid/phasing

RAWDATAFILENAME=B_ribo-seq_manuel_Resting_WT1

fastqc -t 4 -o $OUTPATH/fastqc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz

# TODO create a timestamp folder?

# Trim Adapter
# Only keep relevant RPFs 26-34nt
cutadapt -a AGATCGGAAGAGC -f fastq -e 0.05 -O 3 --quality-base=33 --trim-n -m 20 -M 35 --max-n=0.05 -o $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq <(gzip -dc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz) > $OUTPATH/cutadapt/log/${RAWDATAFILENAME}_trimmed-run-$(date +"%Y-%b-%d_%H-%M-%S").log
 
gzip $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq

fastqc -t 4 -o $OUTPATH/fastqc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz

# Remove Contanminant
bowtie2 --local -k 100 -p 4 --un $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq -x $REFRIBOFILTER -U <(gzip -dc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz) >$OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_filtered.sam 2>$OUTPATH/bowtie2/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log

samtools view -bS $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_filtered.sam >$OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_filtered.bam
rm $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_filtered.sam 

gzip $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq

fastqc -t 4 -o $OUTPATH/fastqc $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz

######## Test Bowtie to remove Contanminant ##########

bowtie -a --best --strata -S --seed 23 -p 4 --chunkmbs 256 --norc --maqerr=60 --un $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_unfiltered.fq $REFRIBOFILTER <(gzip -dc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz) 2>>$OUTPATH/bowtie/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_filtered.sam 

samtools view -bS $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_filtered.sam | samtools sort - -o $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_filtered.bam
samtools index $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_filtered.bam
rm $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_filtered.sam
gzip $OUTPATH/bowtie/${RAWDATAFILENAME}_trimmed_unfiltered.fq

###########################################

# Transcriptome Alignment
bowtie2 --local -p 4 -x $REFTRANSCRIPTOME -U <(gzip -dc $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz) >$OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sam 2>$OUTPATH/bowtie2/log/${RAWDATAFILENAME}_transcriptome-run-$(date +"%Y-%b-%d_%H-%M-%S").log

samtools view -bS $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools idxstats $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sorted.bam.totalreads

samtools view -q 4 -bS $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools idxstats $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >$OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.totalreads
rm $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome.sam 

fastqc -t 4 -o $OUTPATH/fastqc/bowtie2/transcriptome -f bam $OUTPATH/bowtie2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam

# Merge bams for unique mapping and multiple mapping after all samples are processed
samtools merge $OUTPATH/bowtie2/merged_unique.bam $OUTPATH/bowtie2/*_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie2/merged_unique.bam

samtools merge $OUTPATH/bowtie2/merged_multiple.bam $OUTPATH/bowtie2/*_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie2/merged_multiple.bam

######## Test Bowtie2 end-to-end ##########
bowtie2 -p 4 -x $REFTRANSCRIPTOME -U <(gzip -dc $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz) >$OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sam 2>$OUTPATH/bowtie2/end-to-end/log/${RAWDATAFILENAME}_transcriptome-run-$(date +"%Y-%b-%d_%H-%M-%S").log

samtools view -bS $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools idxstats $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sorted.bam.totalreads

samtools view -q 4 -bS $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools idxstats $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >$OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.totalreads
rm $OUTPATH/bowtie2/end-to-end/${RAWDATAFILENAME}_transcriptome.sam 

# Merge bams for unique mapping and multiple mapping after all samples are processed
samtools merge $OUTPATH/bowtie2/end-to-end/merged_unique.bam $OUTPATH/bowtie2/end-to-end/*_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie2/end-to-end/merged_unique.bam

samtools merge $OUTPATH/bowtie2/end-to-end/merged_multiple.bam $OUTPATH/bowtie2/end-to-end/*_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie2/end-to-end/merged_multiple.bam
########################

####### Try Bowtie for transcriptome alignment #######
bowtie -a --best --strata -S -m 100 --seed 23 -p 4 --chunkmbs 256 --norc $REFTRANSCRIPTOME <(gzip -dc $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz) 2>>$OUTPATH/bowtie/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sam

samtools view -bS $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools idxstats $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sorted.bam.totalreads

samtools view -q 4 -bS $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools idxstats $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >$OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.totalreads
rm $OUTPATH/bowtie/${RAWDATAFILENAME}_transcriptome.sam

samtools merge $OUTPATH/bowtie/merged_unique.bam $OUTPATH/bowtie/*_transcriptome_q4.sorted.bam
samtools index $OUTPATH/bowtie/merged_unique.bam

samtools merge $OUTPATH/bowtie/merged_multiple.bam $OUTPATH/bowtie/*_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie/merged_multiple.bam

#######

# Genome Alignment
hisat2 -q -p 4 -x $REFGENOME \
       --known-splicesite-infile $SPLICESITESFILE \
       -U $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz \
       --un-gz $OUTPATH/hisat2/${RAWDATAFILENAME}_nohit.fastq.gz \
       -S $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam 2>$OUTPATH/hisat2/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log

#***********       
#tophat -o $OUTPATH/tophat/test \
#       -g 2 -p 4 --b2-very-sensitive --transcriptome-only --no-novel-juncs \
#       -G $PROTEINCODINGGTF \
#       --transcriptome-index=/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/fasta/protein_coding \
#       $REFGENOME $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz

tophat -p 4 --bowtie1 --no-novel-juncs --GTF $PROTEINCODINGGTF \
       -o $OUTPATH/tophat $REFGENOME $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz 2>$OUTPATH/tophat/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log
#***********

samtools view -bS $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools index $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam
samtools idxstats $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sorted.bam.totalreads

samtools view -q 4 -bS $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam | samtools sort - -o $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools index $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam
samtools idxstats $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats
cut -f3 $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.idxstats | awk '{ sum+=$1} END {print sum}' >$OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam.totalreads
rm $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome.sam 

fastqc -t 4 -o $OUTPATH/fastqc/hisat2/transcriptome -f bam $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam

fastqc -t 4 -o $OUTPATH/fastqc $OUTPATH/hisat2/${RAWDATAFILENAME}_nohit.fastq.gz

# Merge bams for unique mapping and multiple mapping after all samples are processed
samtools merge $OUTPATH/hisat2/merged_unique.bam $OUTPATH/hisat2/*_transcriptome_q4.sorted.bam
samtools index $OUTPATH/hisat2/merged_unique.bam

samtools merge $OUTPATH/hisat2/merged_multiple.bam $OUTPATH/hisat2/*_transcriptome.sorted.bam
samtools index $OUTPATH/hisat2/merged_multiple.bam

### STAR genome ###
# Contanminant Removal
mkdir -p ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/STAR --genomeFastaFiles ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/fasta/ribo-seq_reads-filter.fa --genomeSAindexNbases 5 --genomeChrBinNbits 11

STAR --runThreadN 4 --genomeDir ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/STAR \
     --readFilesIn ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/cutadapt/B_ribo-seq_manuel_Resting_WT1_trimmed.fastq.gz --readFilesCommand gunzip -c \
     --outFileNamePrefix ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/STAR/STAR-contanminant-removal/ \
     --seedSearchLmax 10 \
     --outReadsUnmapped Fastx \
     --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 \
     --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical

# Transcriptome Alignment
mkdir -p ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/STAR --genomeFastaFiles ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/pc_and_lncRNA_transcripts.fa --genomeSAindexNbases 11 --genomeChrBinNbits 12



# Genome Alignment
mkdir -p /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR --genomeFastaFiles ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa --sjdbGTFfile $PROTEINCODINGGTF --sjdbOverhang 34 

#STAR --runThreadN 4 --genomeDir /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR --alignEndsType EndToEnd --readFilesIn $OUTPATH/bowtie2/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz --readFilesCommand gunzip -c --outFilterMismatchNmax 4 --outFilterMultimapNmax 8 --chimScoreSeparation 10 --chimScoreMin 20 --chimSegmentMin 15 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJoverhangMin 500 --outFileNamePrefix "STAR_genome_" --outReadsUnmapped Fastx 
STAR --runThreadN 4 --genomeDir /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR \
     --readFilesIn /Users/huf/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz --readFilesCommand gunzip -c \
     --seedSearchLmax 10 \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/STAR/STAR-genome-test/  \
     --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 4 \
     --outFilterMultimapScoreRange 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax 1 --alignSJoverhangMin 500 

~/BioinformaticsToolShed/STAR/2.5.1b/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR \
--readFilesIn ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/bowtie-contanminant-removal//B_ribo-seq_manuel_Resting_WT1_trimmed_unfiltered.fq.gz --readFilesCommand gunzip -c \
--seedSearchStartLmax 20 \
--sjdbOverhang 34 \
--outReadsUnmapped Fastx \
--outFileNamePrefix ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/STAR/STAR-genome-test10/ \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 \
--outFilterScoreMin 10 --outFilterMatchNmin 15 --outFilterMismatchNmax 2 \
--alignSJoverhangMin 500 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

~/BioinformaticsToolShed/STAR/2.5.1b/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR --readFilesIn ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/bowtie-contanminant-removal//B_ribo-seq_manuel_Resting_WT1_trimmed_unfiltered.fq.gz --readFilesCommand gunzip -c --seedSearchStartLmax 20 --sjdbOverhang 34 --outReadsUnmapped Fastx --outFileNamePrefix ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-09_12-52-44/STAR/STAR-genome-test11/ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outFilterScoreMin 10 --outFilterMatchNmin 15 --outFilterMismatchNmax 2


###################ß

# Transcript selection

# P-site
metagene generate $OUTPATH/plastid/psite/protein_coding \
                  --landmark cds_start \
                  --annotation_files $PROTEINCODINGGTF

psite $OUTPATH/plastid/psite/protein_coding_rois.txt $OUTPATH/plastid/psite/${RAWDATAFILENAME}_hisat2 \
      --min_length 22 --max_length 35 --require_upstream \
      --count_files $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam                 

# Phasing
# TODO how many reads in the given transcripts regions comparing to total reads after alignment
samtools view -hL $PHASINGBED $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >

# Is this counted by offset 12 or 13? Need a check!!!
phase_by_size $OUTPATH/plastid/phasing/${RAWDATAFILENAME}_hisat2 \
              --count_files $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam \
              --annotation_files $PHASINGGTF \
              --annotation_format GTF2 \
              --fiveprime_variable --offset $OUTPATH/plastid/psite/${RAWDATAFILENAME}_hisat2_p_offsets.txt \
              --codon_buffer 5

phase_by_size $OUTPATH/plastid/phasing/test/${RAWDATAFILENAME}_hisat2 \
              --count_files $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam \
              --annotation_files $PHASINGGTF \
              --annotation_format GTF2 \
              --fiveprime_variable --offset $OUTPATH/plastid/psite/${RAWDATAFILENAME}_hisat2_p_offsets.txt \
              --codon_buffer 5 


# ORF calling

#                     
