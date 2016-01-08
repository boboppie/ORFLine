###########################################################
## STEP I - READ PROCESSING AND CONTAMINATION REMOVAL
## 
## a) Clip/trim adaptors
## b) Run bowtie to remove rRNAs
## c) Run bowtie to remove tRNAs
## 
##
## PARAMS to set:
## * Add all executables on PATH environment variable: fastx toolkit, bowtie, tophat, picard
## * Ensure that the variables are set at the top of the script:
## * IGENOMES: path to iGenomes folder containing the Bowtie Indexes for the required species
## * RIBOSOMAL_BOWTIE: Location of rRNA sequences to be searched (in bowtie indexed form)
## * tRNA_BOWTIE: Location of tRNA of sequences to be searched (in bowtie indexed form)
## * PROJECT_PATH: Location of the Project folder within the flowcell (following CASAVA 1.8.2 structure)
## * SAMPLE_NAME: The name of the sample as listed in the sample sheet
###########################################################
IGENOMES="PATH_TO_IGENOMES"
RIBOSOMAL_BOWTIE="PATH_TO_BOWTIE_RIBOSOMAL"
tRNA_BOWTIE="PATH_TO_BOWTIE_tRNA"
PROJECT_PATH="PATH_TO_PROJECT_FOLDER"
SAMPLE_NAME="SAMPLE_NAME"
echo "#Commands used by ARTSeq.sh" > commands.txt
# a) trim/clip adaptors
zcat $PROJECT_PATH/$SAMPLE_NAME/*.fastq.gz | fastx_clipper -a AGATCGGAAGAGCACACGTCT -l 25 -c -n -v -Q33 2>$SAMPLE_NAME.fc.log | fastx_trimmer -Q33 -f 1 2> $SAMPLE_NAME.ft.log > $SAMPLE_NAME.trimmed.fastq
echo "zcat $PROJECT_PATH/$SAMPLE_NAME/*.fastq.gz | fastx_clipper -a AGATCGGAAGAGCACACGTCT -l 25 -c -n -v -Q33 2>$SAMPLE_NAME.fc.log | fastx_trimmer -Q33 -f 1 2> $SAMPLE_NAME.ft.log > $SAMPLE_NAME.trimmed.fastq" > commands.txt

# b) bowtie alignment to remove rRNA contamination
bowtie -p 32 --un=norRNA.Trimmed_$SAMPLE_NAME.fastq $RIBOSOMAL_BOWTIE -l 20 $SAMPLE_NAME.trimmed.fastq 2>> $SAMPLE_NAME.trimmed.stats > $SAMPLE_NAME.rrnaAlignment.fastq
echo "bowtie -p 32 --un=norRNA.Trimmed_$SAMPLE_NAME.fastq $RIBOSOMAL_BOWTIE -l 20 $SAMPLE_NAME.trimmed.fastq 2>> $SAMPLE_NAME.trimmed.stats > $SAMPLE_NAME.rrnaAlignment.fastq" >> commands.txt

# c) bowtie alignment to remove tRNA contamination
echo "bowtie rRNA done"
bowtie -p 32 --un=notRNA_norRNA.Trimmed_$SAMPLE_NAME.nocontam.fastq $tRNA_BOWTIE -l 20 norRNA.Trimmed_$SAMPLE_NAME.fastq  2>> notRNA_norRNA.Trimmed_$SAMPLE_NAME.stats.txt > tRNA.Trimmed_$SAMPLE_NAME.alignment.fastq
echo "bowtie -p 32 --un=notRNA_norRNA.Trimmed_$SAMPLE_NAME.nocontam.fastq $tRNA_BOWTIE -l 20 norRNA.Trimmed_$SAMPLE_NAME.fastq  2>> notRNA_norRNA.Trimmed_$SAMPLE_NAME.stats.txt > tRNA.Trimmed_$SAMPLE_NAME.alignment.fastq" >> commands.txt
echo "bowtie tRNA done"

###########################################################
## STEP 2 - RUN TOPHAT ####
## a) Running tophat
## b) Creating the histogram of trimmed reads
## c) Running picard
## d) Running cufflinks
##
###########################################################

## a) Running tophat
tophat -p 32 --GTF $IGENOMES/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf --no-novel-juncs --output-dir=$SAMPLE_NAME"_tophat"  $IGENOMES/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome notRNA_norRNA.Trimmed_$SAMPLE_NAME.nocontam.fastq 2> $SAMPLE_NAME.tophat.log
echo "tophat -p 32 --GTF $IGENOMES/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf --no-novel-juncs --output-dir=$SAMPLE_NAME"_tophat"  $IGENOMES/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome notRNA_norRNA.Trimmed_$SAMPLE_NAME.nocontam.fastq 2> $SAMPLE_NAME.tophat.log" >> commands.log

echo "tophat done"
# b) Create trimmed read length histograms
samtools view ./$SAMPLE_NAME"_tophat"/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | awk -v v1=0.1 'BEGIN {srand()} !/^$/ { if (rand() <= v1) print length($10)}' | head -n 100000 | sort | uniq -c > ./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME.alignedLengthHistogram.txt
echo "samtools view ./$SAMPLE_NAME"_tophat"/accepted_hits.bam | grep -E '(NM:i:0)|(^@)' | awk -v v1=0.1 'BEGIN {srand()} !/^$/ { if (rand() <= v1) print length($10)}' | head -n 100000 | sort | uniq -c > ./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME.alignedLengthHistogram.txt" >> commands.txt
# run Picard tools for information on distribution of read density across coding,UTR,intronic and intergenic bases
echo "histogram done"
# c) Running picard
java -jar /picard-tools-1.86/CollectRnaSeqMetrics.jar REF_FLAT=$IGENOMES/Homo_sapiens/UCSC/hg19/Annotation/Genes/refFlat.txt.gz STRAND_SPECIFICITY=NONE CHART_OUTPUT=./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME"_metrics.txt" INPUT=./$SAMPLE_NAME"_tophat"/accepted_hits.bam OUTPUT=./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME"_metrics.pdf" 2> $SAMPLE_NAME.picard.log
echo "java -jar /picard-tools-1.86/CollectRnaSeqMetrics.jar REF_FLAT=$IGENOMES/Homo_sapiens/UCSC/hg19/Annotation/Genes/refFlat.txt.gz STRAND_SPECIFICITY=NONE CHART_OUTPUT=./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME"_metrics.txt" INPUT=./$SAMPLE_NAME"_tophat"/accepted_hits.bam OUTPUT=./$SAMPLE_NAME"_tophat"/$SAMPLE_NAME"_metrics.pdf" 2> $SAMPLE_NAME.picard.log" >> commands.txt
echo "picard done"


# d) Running cufflinks
cufflinks --num-threads 16 -G $IGENOMESS/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf $SAMPLE_NAME"_tophat"/accepted_hits.bam 
echo "cufflinks --num-threads 16 -G $IGENOMESS/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf $SAMPLE_NAME"_tophat"/accepted_hits.bam" >> commands.txt
