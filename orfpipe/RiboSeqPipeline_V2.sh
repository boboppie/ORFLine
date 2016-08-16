#!/bin/bash

module load cutadapt/1.9.1
module load fastqc/0.11.4
module load bowtie/1.1.0
module load samtools/1.3

# Parameters
CELLTYPE=B
CONDITION=Resting
AUTHOR=manuel
SAMPLEID=WT1
RAWDATAFILENAME=${CELLTYPE}_ribo-seq_${AUTHOR}_${CONDITION}_${SAMPLEID}

THREADS=4

REFRIBOFILTER=~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/ebwt/ribo-seq_reads-filter
REFRIBOFILTERQUANT=~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/ebwt/ribo-seq_reads-filter-extra-plus-human
REFGENOME=~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome
REFGENOMESTAR=~/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR
REFTRANSCRIPTOME=~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/proteincoding_and_lncRNA-transcript/pc_and_lncRNA_transcripts
SPLICESITESFILE=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.splicesites.txt
PROTEINCODINGGTF=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf
PHASINGGTF=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/plastid.phasing.transcript.gtf
PHASINGBED=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/plastid.phasing.cds.gtf.bed
CHRLIST=~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/chrList.txt


RAWDATAPATH=~/data/orf-discovery/${CELLTYPE}/ribo-seq/${AUTHOR}/${CONDITION}/fastq
OUTPATH=~/out/orf-discovery/${CELLTYPE}/ribo-seq/${AUTHOR}/${CONDITION}

TIMESTAMP=$(date +"%Y-%b-%d_%H-%M-%S")
OUTPATH=$OUTPATH/$TIMESTAMP

mkdir -p $OUTPATH

mkdir -p $OUTPATH/cutadapt/log
mkdir -p $OUTPATH/fastqc/bowtie-contanminant-removal
mkdir -p $OUTPATH/fastqc/bowtie-transcriptome
mkdir -p $OUTPATH/fastqc/star-genome
#mkdir -p $OUTPATH/fastqc/tophat-genome
#mkdir -p $OUTPATH/fastqc/tophat-genome-unmapped
mkdir -p $OUTPATH/bowtie-contanminant-removal/log
mkdir -p $OUTPATH/bowtie-transcriptome/log
mkdir -p $OUTPATH/cufflinks
mkdir -p $OUTPATH/plastid/psite
mkdir -p $OUTPATH/plastid/phasing
mkdir -p $OUTPATH/spectre

###############
fastqc -t $THREADS -o $OUTPATH/fastqc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz

# Trim Adapter
# Only keep relevant RPFs 26-34nt
cutadapt -a AGATCGGAAGAGC -f fastq -e 0.05 -O 3 --quality-base=33 --trim-n -m 25 -M 35 --max-n=0.05 -o $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq <(gzip -dc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz) > $OUTPATH/cutadapt/log/${RAWDATAFILENAME}_trimmed-run-$(date +"%Y-%b-%d_%H-%M-%S").log

gzip $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq

fastqc -t $THREADS -o $OUTPATH/fastqc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz

# Remove Contanminant
bowtie -a --best --strata -S --seed 23 -p $THREADS --chunkmbs 256 --norc --maqerr=60 --un $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq $REFRIBOFILTER <(gzip -dc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz) 2>>$OUTPATH/bowtie-contanminant-removal/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_filtered.sam 

samtools view -bS -F 4 $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_filtered.sam | samtools sort - -o $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_filtered.bam
samtools index $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_filtered.bam
rm $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_filtered.sam

###############################################################
# Contanminant quantification (including extra rRNA annotation)
mkdir -p $OUTPATH/bowtie-contanminant-removal/quantification/log
bowtie -a --best --strata -S --seed 23 -p $THREADS --chunkmbs 256 --norc --maqerr=60 --un $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_unfiltered.fq $REFRIBOFILTERQUANT <(gzip -dc $OUTPATH/cutadapt/${RAWDATAFILENAME}_trimmed.fastq.gz) 2>>$OUTPATH/bowtie-contanminant-removal/quantification/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_filtered.sam 
samtools view -bS -F 4 $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_filtered.sam | samtools sort - -o $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_filtered.bam
samtools index $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_filtered.bam
rm $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_filtered.sam
gzip $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_unfiltered.fq
bowtie -a --best --strata -S --seed 23 -p $THREADS --chunkmbs 256 --norc --un $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_bowtie_transcriptome_unmapped.fq $REFTRANSCRIPTOME <(gzip -dc $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz) 2>>$OUTPATH/bowtie-contanminant-removal/quantification/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_bowtie_transcriptome.sam
gzip $OUTPATH/bowtie-contanminant-removal/quantification/${RAWDATAFILENAME}_bowtie_transcriptome_unmapped.fq

###############################################################

gzip $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq
fastqc -t $THREADS -o $OUTPATH/fastqc/bowtie-contanminant-removal $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz

# Transcriptome Alignment
bowtie -a --best --strata -S --seed 23 -p $THREADS --chunkmbs 256 --norc --un $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_unmapped.fq $REFTRANSCRIPTOME <(gzip -dc $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz) 2>>$OUTPATH/bowtie-transcriptome/log/${RAWDATAFILENAME}-run-$(date +"%Y-%b-%d_%H-%M-%S").log $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sam

gzip $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_unmapped.fq
fastqc -t $THREADS -o $OUTPATH/fastqc/bowtie-transcriptome $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_unmapped.fq.gz

samtools view -bS -F 4 $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sam | samtools sort - -o $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam
samtools index $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam
samtools idxstats $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam >$OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam.totalreads

samtools view -q 255 -b $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sorted.bam >$OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam
samtools index $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam
samtools idxstats $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam >$OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam.idxstats
cut -f3 $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam.totalreads

fastqc -t $THREADS -o $OUTPATH/fastqc/bowtie-transcriptome -f bam $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome_q255.sorted.bam

rm $OUTPATH/bowtie-transcriptome/${RAWDATAFILENAME}_bowtie_transcriptome.sam

# Genome Alignment
#mkdir -p $OUTPATH/tophat-genome/$RAWDATAFILENAME

#tophat -p $THREADS --bowtie1 --no-novel-juncs --GTF $PROTEINCODINGGTF \
#       -o $OUTPATH/tophat/$RAWDATAFILENAME $REFGENOME $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz

#mkdir -p /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR
#STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/STAR --genomeFastaFiles ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa --sjdbGTFfile $PROTEINCODINGGTF --sjdbOverhang 34 

mkdir -p $OUTPATH/star-genome/${RAWDATAFILENAME}
~/BioinformaticsToolShed/STAR/2.5.2a/bin/Linux_x86_64/STAR --runThreadN $THREADS --genomeDir $REFGENOMESTAR \
     --readFilesIn $OUTPATH/bowtie-contanminant-removal/${RAWDATAFILENAME}_trimmed_unfiltered.fq.gz --readFilesCommand gunzip -c \
     --seedSearchStartLmax 20 --sjdbOverhang 34 \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix $OUTPATH/star-genome/${RAWDATAFILENAME}/ \
     --outSAMattributes All --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $THREADS \
     --outFilterScoreMin 10 --outFilterMatchNmin 15 --outFilterMismatchNmax 2 --outFilterMultimapNmax 10 \
     --alignSJoverhangMin 500 --outFilterIntronMotifs RemoveNoncanonicalUnannotated 

# something for small RNA sequencing. REF - http://seqcluster.readthedocs.io/getting_started.html
#STAR --genomeDir $star_index_folder --readFilesIn res/seqs.fastq --alignIntronMax 1  --outFilterMultimapNmax 1000 --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate     

samtools index $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out.bam
samtools view -b -q 255 $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out.bam >$OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out_q255.bam
samtools index $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out_q255.bam

ln -s $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out.bam $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam
ln -s $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out.bam.bai $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam.bai
ln -s $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out_q255.bam $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam
ln -s $OUTPATH/star-genome/$RAWDATAFILENAME/Aligned.sortedByCoord.out_q255.bam.bai $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.bai

samtools idxstats $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam.idxstats
cut -f3 $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome.sorted.bam.totalreads

samtools idxstats $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.idxstats
cut -f3 $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.totalreads

fastqc -t $THREADS -o $OUTPATH/fastqc/star-genome -f bam $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam

# Merge bams for unique mapping and multiple mapping after all samples are processed
samtools merge $OUTPATH/star-genome/merged_q255.bam $OUTPATH/star-genome/*_star_genome_q255.sorted.bam
samtools index $OUTPATH/star-genome/merged_q255.bam

parallel 'samtools view -b $OUTPATH/star-genome/merged_q255.bam {1} >$OUTPATH/star-genome/merged_q255.{1}.bam' :::: $CHRLIST

samtools merge $OUTPATH/star-genome/merged_multiple.bam $OUTPATH/star-genome/*_star_genome.sorted.bam
samtools index $OUTPATH/star-genome/merged_multiple.bam


# cufflinks
cufflinks -p $THREADS -o $OUTPATH/cufflinks -G $PROTEINCODINGGTF $OUTPATH/star-genome/merged_q255.bam

# ##########

# SPECtre is VERY SLOW
python ~/code/github/spectre/SPECtre.py \
        --input $OUTPATH/star-genome/merged_q255.bam \
        --output $OUTPATH/spectre/spectre_output.txt \
        --log $OUTPATH/spectre/spectre_results.log \
        --gtf $PROTEINCODINGGTF \
        --fpkm $OUTPATH/cufflinks/isoforms.fpkm_tracking \
        --floss --orfscore --full
##############################        

# P-site
metagene generate $OUTPATH/plastid/psite/protein_coding \
                  --landmark cds_start \
                  --annotation_files $PROTEINCODINGGTF

psite $OUTPATH/plastid/psite/protein_coding_rois.txt $OUTPATH/plastid/psite/merged_q255_star_genome \
      --min_length 25 --max_length 35 --require_upstream \
      --count_files $OUTPATH/star-genome/merged_q255.bam                

# Phasing
# TODO how many reads in the given transcripts regions comparing to total reads after alignment
# samtools view -hL $PHASINGBED $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam >

# Is this counted by offset 12 or 13? Need a check!!!
phase_by_size $OUTPATH/plastid/phasing/merged_q255_star_genome \
              --count_files $OUTPATH/star-genome/merged_q255.bam \
              --annotation_files $PHASINGGTF \
              --annotation_format GTF2 \
              --fiveprime_variable --offset $OUTPATH/plastid/psite/merged_q255_star_genome_p_offsets.txt \
              --codon_buffer 5

#phase_by_size $OUTPATH/plastid/phasing/test/${RAWDATAFILENAME}_hisat2 \
#              --count_files $OUTPATH/hisat2/${RAWDATAFILENAME}_transcriptome_q4.sorted.bam \
#              --annotation_files $PHASINGGTF \
#              --annotation_format GTF2 \
#              --fiveprime_variable --offset $OUTPATH/plastid/psite/${RAWDATAFILENAME}_hisat2_p_offsets.txt \
#              --codon_buffer 5 


# ORF calling

# Map transcriptome coordinates to genomic coordinates
Rscript ~/code/github/orf-discovery/script/ORFtoBED.R

# Plot ORF length distribution

# Read counts in ORF regions
parallel --dry-run -k --colsep '\t' "samtools view ../../star-genome/merged_q255.chr1.bam {1}:{2}-{3} | wc -l" :::: orfs.bed
bedtools coverage -counts -a orfs.bed -b ../star-genome/merged_q255.bam >orfs_bedtools.counts # this count is wrong!
awk '{if ($13 >0) print}' orfs_bedtools.counts >orfs_bedtools_greaterThanZero.counts
awk '{if ($13 >= 10) print}' orfs_bedtools.counts >orfs_bedtools_greaterThanTen.counts

cut -f1-12 orfs_bedtools_greaterThanZero.counts >orfs_bedtools_greaterThanZero.bed
cut -f11 orfs_bedtools_greaterThanZero.counts | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1' >orfs_bedtools_greaterThanZero.orfwidth
paste -d'\t' orfs_bedtools_greaterThanZero.counts orfs_bedtools_greaterThanZero.orfwidth >orfs_bedtools_greaterThanZero_orfwidth.counts
cut -f1-4,13,14 orfs_bedtools_greaterThanZero_orfwidth.counts >orfs_bedtools_greaterThanZero_region_counts_width
cut -f4 orfs.bed | cut -f1 -d "_" | sort | uniq >txList.txt

awk '{if ($3 == "CDS") print}' ~/ref/GENCODE/mouse/M8/annotation/CHR/comprehensive/annotation.gtf | grep -F -f txList.txt > txList.CDS.gtf
cut -f2 -d";" txList.CDS.gtf | cut -f2 -d"\"" >txList.CDS.gtf.txNames
cut -f1,4,5 txList.CDS.gtf | awk '{print $1"\t"$2-1"\t"$3}' >txList.CDS.bed
paste -d "\t" txList.CDS.bed txList.CDS.gtf.txNames >txList.CDS.txNames.bed

sort-bed txList.CDS.bed >txList.CDS.sorted.bed #461721 records
uniq txList.CDS.sorted.bed >txList.CDS.sorted.uniq.bed #231184 records
sort-bed orfs_bedtools_greaterThanZero_orfwidth.counts >orfs_bedtools_greaterThanZero_orfwidth.counts.sorted
bedops -n orfs_bedtools_greaterThanZero_orfwidth.counts.sorted txList.CDS.sorted.bed >orfs_bedtools_greaterThanZero_orfwidth.counts.sorted.noOverlapWithKnownCDS

sort-bed orfs_bedtools_greaterThanZero.BlocksToCDS.bed >orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed # from https://test.galaxyproject.org Gene BED To Exon/Intron/Codon BED tool

# Remove ORFs with no reads

# for reverse-strand features, counts are reversed relative to genomic coordinates
# REF http://plastid.readthedocs.io/en/latest/examples/count_vector.html?highlight=get_counts
get_count_vectors --annotation_files ~/tmp/plastid/data/orfs_Cxcr4.bed \
                  --annotation_format BED \
                  --count_files ~/tmp/plastid/data/B_resting_chr1.bam \
                  --min_length 25 --max_length 35 \
                  --fiveprime_variable --offset ~/out/orf-discovery/B/ribo-seq/manuel/Resting/2016-May-18_10-44-28/plastid/psite/merged_q255_star_genome_p_offsets.txt \
                  ~/tmp/plastid/data/counts_array

##-----------------------
  502  less txList.CDS.gtf 
  503  less txList.CDS.gtf 
  504  cut -f1,2 -d; txList.CDS.gtf | less
  505  cut -f1,2 -d";" txList.CDS.gtf | less
  506  cut -f1 -d";" txList.CDS.gtf | less
  507  cut -f1-3 -d";" txList.CDS.gtf | less
  508  cut -f2 -d";" txList.CDS.gtf | less
  509  cut -f2 -d";" txList.CDS.gtf | cut -f1 -d""" | less
  510  cut -f2 -d";" txList.CDS.gtf | cut -f1 -d"\"" | less
  511  cut -f2 -d";" txList.CDS.gtf | cut -f2 -d"\"" | less
  512  cut -f2 -d";" txList.CDS.gtf | cut -f2 -d"\"" >txList.CDS.gtf.txNames 
  513  cut -f2 -d";" txList.CDS.gtf | cut -f2 -d"\"" | wc -l
  514  awk '{print $1"\t"$4-1"\t"$5}' txList.CDS.gtf | less
  515  awk '{print $1"\t"$4-1"\t"$5}' txList.CDS.gtf >txList.CDS.bed
  516  awk '{print $1"\t"$4-1"\t"$5}' txList.CDS.gtf | wc -l
  517  lh
  518  paste -d "\t" txList.CDS.bed txList.CDS.gtf.txNames | less
  519  paste -d "\t" txList.CDS.bed txList.CDS.gtf.txNames >txList.CDS.txNames.bed
  520  lh
  521  bed-sort txList.CDS.txNames.bed | less
  522  sort-bed
  523  bedtools
  524  bedops
  525  brew search bedtools
  526  brew search bedops
  527  brew install bedtools bedops
  528  bed-sort txList.CDS.txNames.bed | less
  529  sort-bed txList.CDS.txNames.bed | less
  530  uniq txList.CDS.txNames.bed | sort-bed txList.CDS.txNames.bed | less
  531  uniq txList.CDS.txNames.bed | sort-bed txList.CDS.txNames.bed >txList.CDS.txNames.uniq.sorted.bed 
  532  wc -l txList.CDS.txNames.uniq.sorted.bed 
  533  wc -l txList.CDS.txNames.bed 
  534  uniq txList.CDS.txNames.bed | wc -l
  535  lh
  536  mv txList.CDS.txNames.uniq.sorted.bed txList.CDS.txNames.sorted.bed 
  537  lh
  538  bedops -n 
  539  bedops -n orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | less
  540  bedops -n orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | head -1
  541  head -1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed 
  542  wc -l orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed
  543  bedops -n orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | wc -l
  544  grep ENSMUST00000070533.4 txList.CDS.txNames.sorted.bed 
  545  grep ENSMUST00000070533.4 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed 
  546  bedops -n orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | grep ENSMUST00000070533.4
  547  bedops -e orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | grep ENSMUST00000070533.4
  548  bedops -n 1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | grep ENSMUST00000070533.4
  549  bedops -n 90% orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | grep ENSMUST00000070533.4
  550  bedops -n 90% orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | wc -l
  551  bedops -n 90% orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | head
  552  bedops -n 90% orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed >orfs_tx_overlap_90percent.NOT.bed
  553  lh
  554  cut -f4 orfs_tx_overlap_90percent.NOT.bed | uniq | wc -l
  555  cut -f4 orfs_tx_overlap_90percent.NOT.bed | wc -l
  556  cut -f4 orfs_tx_overlap_90percent.NOT.bed | head
  557  cut -f4 orfs_tx_overlap_90percent.NOT.bed | cut -f1 -d"_" | uniq | wc -l
  558  cut -f4 orfs_tx_overlap_90percent.NOT.bed | cut -f1 -d"_" | uniq |  head
  559  wc -l ENSMUST00000161581ENSMUST00000161581
  560  wc -l txList.CDS.gtf.txNames 
  561  lh
  562  less orfs_tx_overlap_90percent.NOT.bed
  563  grep NSMUST00000027036.10_ORF_1111_1137_0 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed
  564  grep NSMUST00000027036.10 txList.CDS.bed 
  565  grep ENSMUST00000027036.10 txList.CDS.bed 
  566  grep ENSMUST00000027036.10 txList.CDS.gtf
  567  grep ENSMUST00000027036.10 txList.CDS.txNames.bed 
  568  lh
  569  lh
  570  less orfs_tx_overlap_90percent.NOT.bed 
  571  grep ENSMUST00000094273.9_ORF_98_520_1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed
  572  grep ENSMUST00000094273.9 txList.CDS.txNames.bed
  573  grep ENSMUST00000094273.9 orfs_tx_overlap_90percent.NOT.bed 
  574  bedops -n 70% orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed >orfs_tx_overlap_70percent.NOT.bed
  575  lh
  576  wc -l orfs_tx_overlap_70percent.NOT.bed
  577  wc -l orfs_tx_overlap_90percent.NOT.bed
  578  grep ENSMUST00000094273.9 orfs_tx_overlap_70percent.NOT.bed 
  579  bedops -n 1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | less
  580  bedops -n 1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | wc -l
  581  bedops -n 1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed | grep ENSMUST00000094273.9
  582  bedops -n 1 orfs_bedtools_greaterThanZero.BlocksToCDS.sorted.bed txList.CDS.txNames.sorted.bed >orfs_tx_overlap_1bp.NOT.bed                  
      
