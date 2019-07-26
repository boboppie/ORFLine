#!/bin/bash

echo
echo "--------------------------------------------------------------------------------------------"
echo "This script will predict all putative ORFs by scanning transcriptome sequences."
echo
#echo "Prerequisite programmes:"
#echo "parallel"
#echo
echo "Arguments:"
echo "-o|--organism : scientific name, e.g. \"Mus musculus\", default: \"Mus musculus\""
echo "-t|--thread : number of computer cores, default: 4"
echo
echo "Example command: orf-prediction.sh -o \"Mus musculus\" -t 8"
echo "--------------------------------------------------------------------------------------------"
echo

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -o|--organism)
    ORGANISM="$2"
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

if [ "$ORGANISM" == "" ]
then
    ORGANISM="Mus musculus"
fi

if [ "$THREAD" == "" ]
then
    THREAD=4
fi


echo "Input:"
echo "Orgamism - $ORGANISM"
echo "Thread - $THREAD"
echo

REFDIR=./out/ref
GTF=$REFDIR/pipeline_global/annotation.gtf
TRANSCRIPTOME=$REFDIR/pipeline_global/transcriptome.fa


CODEBASE=.
PREDICTIONDIR=./out/prediction

mkdir -p $PREDICTIONDIR

echo "Generating putative ORFs with ATG start codon..."

Rscript $CODEBASE/util/orf-to-bed.R "$ORGANISM" ATG $TRANSCRIPTOME $GTF $PREDICTIONDIR $THREAD

echo "Finished ATG start codon prediction."
echo

echo "Generating putative ORFs with CTG start codon..."

Rscript $CODEBASE/util/orf-to-bed.R "$ORGANISM" CTG $TRANSCRIPTOME $GTF $PREDICTIONDIR $THREAD

echo "Finished CTG start codon prediction."
echo

echo "Generating putative ORFs with GTG start codon..."

Rscript $CODEBASE/util/orf-to-bed.R "$ORGANISM" GTG $TRANSCRIPTOME $GTF $PREDICTIONDIR $THREAD

echo "Finished GTG start codon prediction."
echo

echo "Generating putative ORFs with TTG start codon..."

Rscript $CODEBASE/util/orf-to-bed.R "$ORGANISM" TTG $TRANSCRIPTOME $GTF $PREDICTIONDIR $THREAD

echo "Finished TTG start codon prediction."
echo

echo "Selecting smORFs (<= 100 codon)..."

parallel "cut -f11 $PREDICTIONDIR/orfs_{}.bed | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' >$PREDICTIONDIR/orfs_{}.width" ::: ATG CTG TTG GTG

# 100 codons don't include stop codon, so base length should be 100*3+3=303; what if we extend 100 to 120 codons?
parallel "paste -d'\t' $PREDICTIONDIR/orfs_{}.bed $PREDICTIONDIR/orfs_{}.width | awk '\$13 <= 303' | cut -f1-12 | sort -k4 >$PREDICTIONDIR/orfs_{}.smORFs.bed" ::: ATG CTG TTG GTG

#nice -5 parallel "paste -d'\t' orfs_{}.bed orfs_{}.width | awk '\$13 <= 363' | cut -f1-12 | sort -k4 >orfs_{}.smORFs_120.bed" ::: ATG CTG TTG GTG

echo "Finished ORF prediction."
echo
