#!/bin/bash

echo
echo "----------------------------------------------------------------------------------------"
echo "This script will process the raw RNA-Seq fastq file."
echo
echo "Steps:"
echo "1. QC"
echo "2. Adapter/quality Trimming"
echo "3. Reference genome alignment"
echo "4. Transcript expression quantification"
echo
echo "Prerequisite programmes:"
echo "FastQC"
echo "Trim Galore"
echo "Samtools"
echo "STAR"
echo "StringTie"
echo
echo "Arguments:"
echo "-f|--fastq : e.g. path to a FSATQ file. If paried-end, expects"
echo "             files to be supplied in a pairwise fashion, e.g. file_1.fq.gz file_2.fq.gz"
echo "-p|--paired : if the reads are pairec-end, either true or false, default: false"
echo "-t|--thread : number of computer cores, default: 4"
echo
echo "Example command: rnaseq-process.sh -f ./out/data/rna-seq/rna.fastq.gz -t 4"
echo "----------------------------------------------------------------------------------------"
echo

echo
echo "----------------------------------------------------------------------------------------"
echo "Checking the availability of tools needed for the pipeline..."
echo "----------------------------------------------------------------------------------------"
echo

# Declare a string array with type
declare -a CmdArray=("samtools" "fastqc" "trim_galore" "stringtie" "STAR")
 
# Read the array values with space
for cmd in "${CmdArray[@]}"; do
  if ! [ -x "$(command -v $cmd)" ]; then
  	echo
    echo "Error: $cmd is not installed." >&2
    exit 1
  else
    echo "$cmd is available."  
  fi
done


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--fastq)
    FASTQ="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--paired)
    PAIRED="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--thread)
    THREAD="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

DATADIR=./out/data

if [ "$FASTQ" == "" ]
then
    FASTQ=$DATADIR/rna-seq/test.fastq.gz
fi

if [ "$PAIRED" == "" ]
then
    PAIRED=false
fi

if [ "$THREAD" == "" ]
then
    THREAD=4
fi

echo
echo "Input:"
echo "FASTQ - $FASTQ"
echo "Paired-end - $PAIRED"
echo "Thread - $THREAD"
echo

t=${FASTQ%.fastq.gz}
NAME=${t##*/}

REFDIR=./out/ref
REFGENOMESTAR=$REFDIR/pipeline_global/star/100bp
GTF=$REFDIR/pipeline_global/annotation.gtf

OUTPATH=./out/process/rna

mkdir -p $OUTPATH/trim_galore


echo "Running Trim Galore with FastQC..."

if [ "$PAIRED" == "true" ]
then
    trim_galore -q 33 --fastqc --trim-n -e 0.1 --stringency 3 --paired  $FASTQ -o $OUTPATH/trim_galore
else 
    trim_galore -q 33 --fastqc --trim-n -e 0.1 --stringency 3 $FASTQ -o $OUTPATH/trim_galore
fi

echo "Finished Trim Galore..."
echo


echo "Running STAR for genome alignment..."

SAM_ATTR="NH HI NM MD AS"
MISMATCH_MAX=1

mkdir -p $OUTPATH/star-genome/$NAME

STAR --runThreadN $THREAD --genomeDir $REFGENOMESTAR \
         --readFilesIn $OUTPATH/trim_galore/${NAME}*trimmed.fq.gz --readFilesCommand zcat \
         --alignEndsType EndToEnd \
         --outFilterMismatchNmax $MISMATCH_MAX \
         --outReadsUnmapped Fastx --outFileNamePrefix $OUTPATH/star-genome/$NAME/ \
         --outSAMattributes $SAM_ATTR --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $THREAD  

#if [ "$PAIRED" == "true" ]
#then
#    STAR --runThreadN $THREADS --genomeDir $REFGENOMESTAR100SE \
#         --readFilesIn $OUTPATH/trim_galore/${RAWDATAFILENAME}_R1_trimmed.fq.gz $OUTPATH/trim_galore/${RAWDATAFILENAME}_R2_trimmed.fq.gz --readFilesCommand zcat \
#         --alignEndsType EndToEnd \
#         --outFilterMismatchNmax $MISMATCH_MAX \
#         --outReadsUnmapped Fastx --outFileNamePrefix $OUTPATH/star/${RAWDATAFILENAME}/ \
#         --outSAMattributes $SAM_ATTR --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $THREADS
#else 
     
#fi


# Note: align to transcriptome, add --quantMode TranscriptomeSAM 

gzip $OUTPATH/star-genome/$NAME/Unmapped.out.mate* &

samtools index $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam &
samtools view -b -q 255 $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam >$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam
samtools index $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam

#fastqc -t $NCORE -o $OUTPATH/fastqc/star-genome -f bam $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam &

ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam.bai $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam.bai
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam.bai $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam.bai


echo "Finished STAR."
echo


echo "Running StringTie for gene/transcript expression estimation..."

mkdir -p $OUTPATH/stringtie/$NAME

stringtie $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam  -p $THREAD -G $GTF -eB -o $OUTPATH/stringtie/$NAME/expressed.gtf -A $OUTPATH/stringtie/$NAME/gene_abund.tab -C $OUTPATH/stringtie/$NAME/cov_refs.gtf

cat $OUTPATH/stringtie/$NAME/expressed.gtf | awk '{if($3=="transcript") print}' | cut -f9 | cut -f4,6,10,12 -d" " | sed 's/"//g' | sed 's/;//g' | tr ' ' \\t | sort >$OUTPATH/stringtie/expressed.tsv

# FPKM cutoff set to 0.5

cat $OUTPATH/stringtie/expressed.tsv | awk '$3 > 0.5' | cut -f1 | sort | uniq >$OUTPATH/stringtie/active_transcript_Ids.txt 

# StringTie wrongly estimate expression for transcripts, so if a transcript pass the expression cutoff, take all transcript from the same gene
LL_ALL=C fgrep -f $OUTPATH/stringtie/active_transcript_Ids.txt $OUTPATH/stringtie/expressed.tsv | cut -f2 | sort | uniq >$OUTPATH/stringtie/active_transcript_geneName.txt
LL_ALL=C fgrep -f $OUTPATH/stringtie/active_transcript_geneName.txt $OUTPATH/stringtie/expressed.tsv | cut -f1 | sort | uniq >$OUTPATH/stringtie/active_transcript_gene_transcriptAll.txt

mv $OUTPATH/stringtie/active_transcript_Ids.txt $OUTPATH/stringtie/active_transcript_Ids.txt_
ln -s ../../../../$OUTPATH/stringtie/active_transcript_gene_transcriptAll.txt $OUTPATH/stringtie/active_transcript_Ids.txt

echo "Finished StringTie."
echo

echo "Finished RNA-Seq processing."
echo
