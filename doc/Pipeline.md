## Data
~/data/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastq

## Annotation
* Whole genome
~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa

* tRNA, rRNA, Mt_rRNA, Mt_tRNA, snRNA, snoRNA, misc_RNA, miRNA
/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter.fa

## Workflow
* Software
OS X        v10.11.2
samtools    v1.2        (brew)
FastQC      v0.11.4     (direct download for linux)
cutadapt    v1.9.1      (pip)
trim_galore v0.4.1      (direct download)
bowtie      v1.1.2      (direct download pre-compiled for macosx, brew installation produces rubbish)
bowtie2     v2.2.6      (brew)
TopHat      v2.0.14     (brew)
cufflinks   v2.2.1      (direct download pre-compiled)
kallisto    v0.42.4     (brew)
plastid     v0.4.4      (pip)
bedops      v2.4.15     (brew)
hisat2      v2.0.1-beta (direct download pre-compiled)

1. Trim adaptor (QC), Sanger format (Phred+33)
###trim_galore -q 15 --phred33 --stringency 3 -e 0.05 --fastqc --illumina -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/trim_galore ~/data/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastq/BMDM_ribo-seq_manuel_LPS_TTP_WT1.fastq.gz 

cutadapt -a AGATCGGAAGAGC -f fastq -e 0.05 -O 3 --quality-base=33 --trim-n -m 20 -M 35 --max-n=0.05 -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.fastq <(gzip -dc ~/data/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastq/BMDM_ribo-seq_manuel_LPS_TTP_WT1.fastq.gz) > ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.report.txt
 
gzip ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.fastq

fastqc -t 4 -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastqc ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.fastq.gz

2. Map reads to ribo-seq_reads-filter.fa by bowtie
###bowtie --verbose -a --best --strata -S -m 100 --seed 23 -p 4 --chunkmbs 256 --un ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq  ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter <(gzip -dc ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.fastq.gz) | samtools view -bS - >~/out/orf-discovery/BMDtr '\n' ' 'M/ribo-seq/manuel/LPS/bowtie/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_filtered.bam

bowtie2 --local -k 100 --seed 23 -p 4 --un ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq -x ~/ref/GENCODE/mouse/M8/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter -U <(gzip -dc ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cutadapt/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed.fastq.gz) | samtools view -bS - >~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_filtered.bam
gzip ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq

3. Map unfiltered reads to genome/transcriptome by  , detect novel transcripts by Cufflinks
# DO NOT STREAM GTF
# Align to genome
tophat -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome -g 64 -p 4 \ 
       -G ~/ref/GENCODE/mouse/M8/annotation/CHR/filtered/filtered.gtf \ 
       --b2-very-sensitive ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome \ 
       ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq.gz

#to resume:
#tophat -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome -g 64 -p 4 -R ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome -G ~/ref/GENCODE/mouse/M8/annotation/CHR/filtered/filtered.gtf --b2-very-sensitive ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq.gz

samtools index ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome/accepted_hits.bam
samtools view -q 10 -b ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome/accepted_hits.bam > ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome/accepted_hits.q10.bam # mapping quality filtering by lh3 on biostar
samtools index ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome/accepted_hits.q10.bam
fastqc -t 4 -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastqc/tophat/genome -f bam ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/genome/accepted_hits.q10.bam

# Align to transcriptome (cDNA, introns excluded)

################ DO NOT RUN
# Building transcriptome data files (NOTE: tophat only takes relavent path???)
# cd /ref/GENCODE/mouse/M8/annotation/CHR/protein_coding
# tophat -G ~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.gtf --transcriptome-index=tophat-index/protein_coding ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome
# cut -f2 -d" " tophat-index/protein_coding.fa.tlst | sort | uniq >protein_coding.fa.tlst.ids
# diff -y tophat-index/protein_coding.fa.tlst.ids ../../../fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.tlst.ids 
# cd ~
################

# Create a gtf based on fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.tlst.ids, see download.sh
# parallel 'LANG=C grep -F {} ~/ref/GENCODE/mouse/M8/annotation/CHRcomprehensive/annotation.gtf >~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.gtf.{}' :::: ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.tlst.ids
# cat ~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.gtf.* >~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.gtf
# rm ~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.gtf.*

tophat -o /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome \ 
       -g 2 -p 4 --b2-very-sensitive --transcriptome-only --no-novel-juncs \ 
       -G /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf \ 
       --transcriptome-index=/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/fasta/protein_coding \ 
       /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome \ 
       /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq.gz

samtools index ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.bam
samtools view -q 4 -b ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.bam > ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam # mapping quality filtering by lh3 on biostar, however as Siomon suggested, 4 and above should be uniq primary mapping
samtools index ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam
fastqc -t 4 -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastqc/tophat/transcriptome -f bam ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam

# From Simon - mapping reads to genome (HISAT2) and filter out low mapping quality reads rather than using tophat transcriptome only mode
extract_splice_sites.py /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf  > ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/splicesites.txt
hisat2-build /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome
hisat2 -q -p 4 -x /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome --known-splicesite-infile ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/splicesites.txt --seed 23 -U /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq.gz -S /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome.sam

samtools view -q 4 -bS /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome.sam > /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.bam # 67049709 -> 23869793
samtools sort /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.bam -o /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam
samtools index /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam
fastqc -t 4 -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/fastqc/hisat2/transcriptome -f bam /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam

# to calculate read length mean and sd, using samtools stats and a R func - http://stackoverflow.com/questions/22644481/r-computing-mean-median-variance-from-file-with-frequency-distribution
# samtools stats ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/accepted_hits.q10.bam
#      mean    median       var        sd 
# 29.207389 29.000000  3.817538  1.953852 

######################### WARNING - Applying Cufflinks and kallisto to Ribo-seq data may break the underline assumeption of the read distribution, also Cufflinks is not realiable on novel transcript construncton

cufflinks -o ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks -p 4 --seed 23 -g ~/ref/GENCODE/mouse/M8/annotation/CHR/filtered/filtered.gtf ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/accepted_hits.q10.bam

# filtered out transcripts that have coverage (Estimate for the absolute depth of read coverage across the tranascript) < 0.4
awk '{if ($9 >= 0.4) print}' ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/isoforms.fpkm_tracking >~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/isoforms.fpkm_tracking.coverageabovepoint4

# Quantify abundances of transcripts by kallisto (or Sailfish, Salmon, etc.)
gtf_juncs ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/transcripts.gtf  > ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts.juncs

gtf_to_fasta --min-anchor 8 --splice-mismatches 0 --min-report-intron 50 --max-report-intron 500000 --min-isoform-fraction 0.15 --output-dir /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/ --max-multihits 64 --max-seg-multihits 128 --segment-length 25 --segment-mismatches 2 --min-closure-exon 100 --min-closure-intron 50 --max-closure-intron 5000 --min-coverage-intron 50 --max-coverage-intron 20000 --min-segment-intron 50 --max-segment-intron 500000 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --read-realign-edit-dist 3 --max-insertion-length 3 --max-deletion-length 3 -z gzip -p4 --gtf-annotations ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/transcripts.gtf --gtf-juncs ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts.juncs --no-closure-search --no-coverage-search --no-microexon-search ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/transcripts.gtf  /Users/huf/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts.fa > /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/g2f.log

awk '/^>/{print ">"$2; next}{print}' /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts.fa >/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts_renamed.fa

kallisto index -i /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts_renamed.idx /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts_renamed.fa

kallisto quant -i /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/cufflinks_transcripts_renamed.idx -o /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto -b 100 --seed=23 --single -l 200 -s 80 -t 4 ~/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/bowtie2/BMDM_ribo-seq_manuel_LPS_TTP_WT1_trimmed_unfiltered.fq.gz

# awk '{if($4>0 && $4<1)print}' abundance.tsv > abundance_low.tsv
# awk '{if($4>1)print}' abundance.tsv > abundance_high.tsv
# join -j 1 <(sort -k1 abundance_high.tsv) <(sort -k1 ../cufflinks/isoforms.fpkm_tracking.belowpoint4) | cut -f1,4,5,11,13,14 -d" " > kallisto_high_cufflinks_low_overlap.txt

# filtered out transcripts that have est_counts < 5
awk '{if($4>=5)print}' /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/abundance.tsv > /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/abundance_estcountsabove5.tsv

# Take union of cufflinks and kallisto transcript ids
cat /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/abundance_estcountsabove5.tsv /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/isoforms.fpkm_tracking.abovepoint4 | cut -f1 | egrep 'ENSMUS|CUFF' | sort | uniq >/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/kallisto_cufflinks_transcripts_id.txt

####################################################

# and gtf
grep -f /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/kallisto_cufflinks_transcripts_id.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/cufflinks/transcripts.gtf >/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/kallisto/kallisto_cufflinks_transcripts.gtf

4. P-site/metagene profile/triplet periodicity (QC)

mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite
mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript
mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/phasing
mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/expression

# Use transcriptome bam
# Piste
metagene generate /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/protein_coding \ 
                  --landmark cds_start \ 
                  --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf

mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/hisat2
mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/tophat
psite /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/protein_coding_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/tophat/riboprofile_tophat_transcriptome --min_length 22 --max_length 33 --require_upstream --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam

psite /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/protein_coding_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/hisat2/riboprofile_tophat_transcriptome --min_length 22 --max_length 32 --require_upstream --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam

# Metagene analysis
metagene generate /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/protein_coding  --landmark cds_start \
                  --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf \
                  --downstream 200

metagene count /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/protein_coding_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome \
               --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
               --fiveprime --offset 12

metagene chart /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome \
               /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome_metagene_profile.txt \
               --landmark "start codon" \
               --title "Metagene demo"

mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2
metagene generate /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2/protein_coding  --landmark cds_start \
                  --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf \
                  --downstream 30

metagene count /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2/protein_coding_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2/riboprofile_hisat2_transcriptome \
               --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam \
               --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/hisat2/riboprofile_tophat_transcriptome_p_offsets.txt 

metagene chart /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2/riboprofile_hisat2_transcriptome \
               /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/hisat2/riboprofile_hisat2_transcriptome_metagene_profile.txt \
               --landmark "start codon" \
               --title "Metagene profile"              

## Or give an offsets file
metagene count /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/protein_coding_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome_fiveprime_variable \
               --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
               --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/riboprofile_tophat_transcriptome_p_offsets.txt

metagene chart /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome_fiveprime_variable \
               /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/riboprofile_tophat_transcriptome_fiveprime_variable_metagene_profile.txt \
               --landmark "start codon" \
               --title "Metagene demo"

### Try a single transcript:
ENSID=ENSMUST00000110082
grep $ENSID ~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf >~/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/${ENSID}.gtf

metagene generate /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID} --landmark cds_start \
                  --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/${ENSID}.gtf \
                  --downstream 200

metagene count /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID}_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID} \
               --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
               --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/riboprofile_tophat_transcriptome_p_offsets.txt

metagene chart /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID} \
               /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID}_metagene_profile.txt \
               --landmark "start codon" \
               --title "${ENSID}"

mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/hisat2
metagene count /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/${ENSID}_rois.txt /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/hisat2/${ENSID} \
               --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam \
               --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/hisat2/riboprofile_tophat_transcriptome_p_offsets.txt 

metagene chart /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/hisat2/${ENSID} \
               /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/metagene/transcript/hisat2/${ENSID}_metagene_profile.txt \
               --landmark "start codon" \
               --title "${ENSID}"
###

# Phasing
awk '{if($3=="exon" || $3=="CDS" || $3=="start_codon" || $3=="stop_codon"){print $0}}' protein_coding.filtered.noSelenocysteine.gtf >protein_coding.filtered.noSelenocysteine.plastid.gtf
awk '{if($0~"_level \"(1|2|3)\";"){print $0}}' protein_coding.filtered.noSelenocysteine.plastid.gtf >protein_coding.filtered.noSelenocysteine.plastid.supportLevel1To3.gtf
awk '{if($0!~"tag \"(cds_end_NF|cds_start_NF)\";"){print $0}}' protein_coding.filtered.noSelenocysteine.plastid.gtf >protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.gtf
grep -v non_stop_decay protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.gtf >protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.gtf
awk '{if($0!~"_level \"NA\";"){print $0}}' protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.gtf >protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.gtf
awk '{if($0~"level (1|2);"){print $0}}' protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.gtf >protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.level12.gtf
grep -v Ntan1-004 protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.level12.gtf >protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.level12.noNtan1-004.gtf

phase_by_size /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/phasing/riboprofile_tophat_transcriptome \
                --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
                --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf\
                --annotation_format GTF2 \
                --fiveprime --offset 12 \
                --codon_buffer 5 \
                --min_length 28 --max_length 30

phase_by_size /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/phasing/riboprofile_tophat_transcriptome \
                --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam \
                --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.plastid.tagFiltered.nonStopDecay.noSupportLevelNA.level12.noNtan1-004.gtf \
                --annotation_format GTF2 \
                --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/hisat2/riboprofile_tophat_transcriptome_p_offsets.txt \
                --codon_buffer 5  

# Expression
counts_in_region /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/expression/riboprofile_tophat_transcriptome_counts.txt --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
                 --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf \
                 --fiveprime --offset 12

counts_in_region /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/expression/riboprofile_tophat_transcriptome_counts_fiveprime_variable.txt --count_files /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/tophat/transcriptome/accepted_hits.q10.bam \
                 --annotation_files /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf \
                 --fiveprime_variable --offset /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/plastid/psite/riboprofile_tophat_transcriptome_p_offsets.txt                


5. Transcripts selection

# ORF-RATER
pip install git+https://github.com/PyTables/PyTables

mkdir -p /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/ORF-RATER/

# gtf2bed < protein_coding.filtered.noSelenocysteine.gtf > protein_coding.filtered.noSelenocysteine.bed # not parsing correctly
# awk '{if ($8 == "transcript") print}' protein_coding.filtered.noSelenocysteine.bed >protein_coding.filtered.noSelenocysteine.transcript.bed

fgrep -w transcript /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.gtf | sed 's/[";]//g;' | awk '{OFS="\t"; print $1, $4-1,$5,$12,0,$7,$18,$14,$10}' >/Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.transcript.bed # REF http://onetipperday.sterding.com/2012/08/convert-bed-to-gtf.html

python ~/code/github/ORF-RATER/prune_transcripts.py --inbed /Users/huf/ref/GENCODE/mouse/M8/annotation/CHR/protein_coding/protein_coding.filtered.noSelenocysteine.transcript.bed --outbed /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/ORF-RATER/transcripts.bed --summarytable /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/ORF-RATER/summary.txt --minlen 28 --maxlen 30 -vv -p 2 -f --keeptempfiles ~/ref/GENCODE/mouse/M8/fasta/genome/CHR/GRCm38.CHR.genome.fa /Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/hisat2/transcriptome_q4.sorted.bam >/Users/huf/out/orf-discovery/BMDM/ribo-seq/manuel/LPS/ORF-RATER/log/run-$(date +%s).log

Upstream open reading frames (uORFs) are prevalent in the human transcriptome and may negatively regulate the abundance of canonically encoding proteins through the promotion of mRNA decay and competitive expression, among other mechanisms. uORFs are conserved across species and have been annotated to genes with diverse biological functions, including but not limited to oncogenes, cell cycle control and differentiation, and stress response. As such, the aberrant expression of certain uORFs has been implicated in the development and progression of various diseases. Therefore, the positive identification and validation of uORFs as translational products is critical for understanding their role in complex biological processes and disease etiology. Where mRNA-Seq has been used to approximate the transcriptomic content of a cell, or group of cells, the recently developed method of sequencing ribosome-protected fragments aims to profile the translational landscape of a sample. In concert, various algorithms have been developed to differentiate coding transcripts from non-coding transcripts based on the alignment of ribosome-protected fragments to a reference transcriptome. We have developed a classification algorithm based on the magnitude of coherence between the aligned ribosome profiling reads and tri-nucleotide periodic signal inherent to protein-coding sequences.

