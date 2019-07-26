#!/bin/bash

echo
echo "--------------------------------------------------------------------------------------------"
echo "This script will process the raw Ribo-Seq fastq file."
echo
echo "Steps:"
echo "1. QC"
echo "2. Adapter/quality trimming"
echo "3. Contanminant removal"
echo "4. Reference genome alignment"
echo "5. P-site calling"
echo
echo "Prerequisite programmes:"
echo "gzip"
echo "FastQC"
echo "Trim Galore"
echo "Bowtie"
echo "Samtools"
echo "STAR"
echo "plastid"
echo
echo "Arguments:"
echo "-f|--fastq : e.g. path to a FSATQ file"
echo "-a|--adapter : adapter sequence, default: Trim Galore auto-detection"
echo "-t|--thread : number of computer cores, default: 4"
echo
echo "Example command: riboseq-process.sh -f ./out/data/ribo-seq/ribo.fastq.gz -a AAAAAAAAAAA -t 4"
echo "--------------------------------------------------------------------------------------------"
echo

echo
echo "--------------------------------------------------------------------------------------------"
echo "Checking the availability of tools needed for the pipeline..."
echo "--------------------------------------------------------------------------------------------"
echo

# Declare a string array with type
declare -a CmdArray=("samtools" "bowtie" "fastqc" "trim_galore" "psite" "STAR" "gzip")
 
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
    -a|--adapter)
    ADAPTER="$2"
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
    FASTQ=$DATADIR/ribo-seq/test.fastq.gz
fi

if [ "$THREAD" == "" ]
then
    THREAD=4
fi

echo
echo "Input:"
echo "FASTQ - $FASTQ"
echo "Adapter - $ADAPTER"
echo "Thread - $THREAD"
echo

t=${FASTQ%.fastq.gz}
NAME=${t##*/}

REFDIR=./out/ref
REFRIBOFILTER=$REFDIR/pipeline_global/contaimination/bowtie_index
REFGENOMESTAR=$REFDIR/pipeline_global/star/50bp
PROTEINCODINGGTF=$REFDIR/pipeline_global/pc.gtf.gz

OUTPATH=./out/process/ribo

mkdir -p $OUTPATH/trim_galore
mkdir -p $OUTPATH/bowtie-contanminant-removal
mkdir -p $OUTPATH/plastid/psite
mkdir -p $OUTPATH/plastid/phasing


echo "Running Trim Galore with FastQC..."

if [ "$ADAPTER" == "" ]
then
   trim_galore -j $THREAD -q 33 --fastqc --trim-n -e 0.1 --stringency 3 $FASTQ -o $OUTPATH/trim_galore # max and min read length: --length 25 --max_length 35
else 
   trim_galore -j $THREAD -q 33 -a $ADAPTER --fastqc --trim-n -e 0.1 --stringency 3 $FASTQ -o $OUTPATH/trim_galore
fi

echo "Finished Trim Galore..."
echo


# Remove Contanminant
echo "Running Bowtie to remove contanminant..."

bowtie -a --best --strata -S --seed 23 -p $THREAD --chunkmbs 256 --norc --maqerr=60 --un $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_unfiltered.fq $REFRIBOFILTER <(gzip -dc $OUTPATH/trim_galore/${NAME}_trimmed.fq.gz) $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_filtered.sam 

samtools view -bS -F 4 $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_filtered.sam | samtools sort - -o $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_filtered.bam
samtools index $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_filtered.bam
rm $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_filtered.sam

gzip $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_unfiltered.fq

echo "Finished Bowtie."
echo


echo "Running STAR for genome alignment..."

SAM_ATTR="NH HI NM MD AS"
ALIGNINTRON_MIN=20
ALIGNINTRON_MAX=10000
MISMATCH_MAX=1
MISMATCH_NOVERL_MAX=0.04
FILTER_TYPE=BySJout

mkdir -p $OUTPATH/star-genome/$NAME


STAR --runThreadN $THREAD --genomeDir $REFGENOMESTAR \
     --readFilesIn $OUTPATH/bowtie-contanminant-removal/${NAME}_trimmed_unfiltered.fq.gz --readFilesCommand zcat \
     --outReadsUnmapped Fastx --outFileNamePrefix $OUTPATH/star-genome/$NAME/ \
     --alignIntronMin $ALIGNINTRON_MIN --alignIntronMax $ALIGNINTRON_MAX --alignEndsType EndToEnd \
     --outFilterMismatchNmax $MISMATCH_MAX --outFilterMismatchNoverLmax $MISMATCH_NOVERL_MAX \
     --outFilterType $FILTER_TYPE --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
     --outSAMattributes $SAM_ATTR --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $THREAD

# Note: align to transcriptome, add --quantMode TranscriptomeSAM 

gzip $OUTPATH/star-genome/$NAME/Unmapped.out.mate1 &

samtools index $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam &
samtools view -b -q 255 $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam >$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam
samtools index $OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam

#fastqc -t $THREAD -o $OUTPATH/fastqc/star-genome -f bam $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam &

ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out.bam.bai $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam.bai
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam
ln -s ../../../../$OUTPATH/star-genome/$NAME/Aligned.sortedByCoord.out_q255.bam.bai $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam.bai

#(samtools idxstats $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam >$OUTPATH/star-genome/${NAME}_star_genome.sorted.bam.idxstats; cut -f3 $OUTPATH/star-genome/${NAME}_star_genome.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/star-genome/${NAME}_star_genome.sorted.bam.totalreads) &

#samtools idxstats $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.idxstats
#cut -f3 $OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.idxstats | paste -sd+ - | bc >$OUTPATH/star-genome/${RAWDATAFILENAME}_star_genome_q255.sorted.bam.totalreads

#----------------------

# Merge bams for unique mapping and multiple mapping after all samples are processed
#(samtools merge $OUTPATH/star-genome/merged_q255.bam $OUTPATH/star-genome/*_star_genome_q255.sorted.bam; samtools index $OUTPATH/star-genome/merged_q255.bam) &

#parallel 'samtools view -b $OUTPATH/star-genome/merged_q255.bam {1} >$OUTPATH/star-genome/merged_q255.{1}.bam' :::: $CHRLIST

#samtools merge $OUTPATH/star-genome/merged_multiple_withUnmapped.bam $OUTPATH/star-genome/*_star_genome.sorted.bam
#samtools view -b -F 4 $OUTPATH/star-genome/merged_multiple_withUnmapped.bam >$OUTPATH/star-genome/merged_multiple.bam
#samtools index $OUTPATH/star-genome/merged_multiple.bam
#rm $OUTPATH/star-genome/merged_multiple_withUnmapped.bam

#samtools view -h merged_q255.bam | awk 'length($10) >= 28 && length($10) <= 30 || $1 ~ /^@/' | samtools view -bS - > merged_q255_28-30nt.bam3
#samtools index merged_q255_28-30nt.bam

echo "Finished STAR."
echo

   
echo "Running metagene analysis..."

metagene generate $OUTPATH/plastid/psite/protein_coding \
                  --landmark cds_start \
                  --annotation_files $PROTEINCODINGGTF

echo "Finished metagene analysis."
echo

echo "Running psite calling..."                  

psite $OUTPATH/plastid/psite/protein_coding_rois.txt $OUTPATH/plastid/psite/${NAME}_q255_star_genome \
      --min_length 25 --max_length 35 --require_upstream \
      --count_files $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam

echo "Finished psite calling."
echo


echo "Running phase calling..."

phase_by_size $OUTPATH/plastid/psite/protein_coding_rois.txt $OUTPATH/plastid/phasing/${NAME}_q255_star_genome \
              --count_files $OUTPATH/star-genome/${NAME}_star_genome_q255.sorted.bam \
              --fiveprime_variable --offset $OUTPATH/plastid/psite/${NAME}_q255_star_genome_p_offsets.txt \
              --codon_buffer 5

echo "Finished phase calling."
echo

echo "Finished Ribo-Seq processing."
echo
