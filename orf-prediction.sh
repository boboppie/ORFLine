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
echo "-c|--codon : a string of start codons, default: \"ATG CTG GTG TTG\""
echo "-u|--upper : upper limit of smORF length (codons), default: 100"
echo
echo "Example command: orf-prediction.sh -o \"Mus musculus\" -t 8 -c \"ATG CTG GTG TTG\" -u 100"
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
    -c|--codon)
    CODONLIST="$2"
    shift # past argument
    shift # past value
    ;;
    -u|--upper)  
    UPPER="$2"
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

if [ "$CODONLIST" == "" ]
then
    CODONLIST="ATG CTG GTG TTG"
fi

if [ "$UPPER" == "" ]
then
    UPPER=100
fi


echo "Input:"
echo "Orgamism - $ORGANISM"
echo "Thread - $THREAD"
echo "Start codons - $CODONLIST"
echo "smORF length upper limit - $UPPER"
echo

REFDIR=./out/ref
GTF=$REFDIR/pipeline_global/annotation.gtf
TRANSCRIPTOME=$REFDIR/pipeline_global/transcriptome.fa


CODEBASE=.
PREDICTIONDIR=./out/prediction

mkdir -p $PREDICTIONDIR

parallel -j1 "echo \"Generating putative ORFs with {} start codon...\"; Rscript $CODEBASE/util/orf-to-bed.R \"$ORGANISM\" {} $TRANSCRIPTOME $GTF $PREDICTIONDIR $THREAD; echo \"Finished {} start codon prediction.\"; echo" ::: $CODONLIST

echo "Selecting smORFs (<= $UPPER codons)..."

parallel "cut -f11 $PREDICTIONDIR/orfs_{}.bed | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' >$PREDICTIONDIR/orfs_{}.width" ::: $CODONLIST

# 100 codons don't include stop codon, so base length should be 100*3+3=303
UPPERINNT=$((UPPER*3+3))

parallel "paste -d'\t' $PREDICTIONDIR/orfs_{}.bed $PREDICTIONDIR/orfs_{}.width | awk '\$13 <= $UPPERINNT' | cut -f1-12 | sort -k4 >$PREDICTIONDIR/orfs_{}.smORFs.bed" ::: $CODONLIST

echo "Finished ORF prediction."
echo
