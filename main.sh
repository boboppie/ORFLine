#!/bin/bash

# Master script

# Parse argument, REF - https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
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
    -r|--release)
    RELEASE="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--file)
    RELEASE="$2"
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

LOGDIR=./out/log
mkdir -p $LOGDIR

echo 'You can also find output messages in ./out/log'
echo

# Step 1 - check dependencies
bash ./module-check.sh 2>&1 | tee $LOGDIR/module-check.log

# Step 2 - download sample data

bash ./data-download.sh 2>&1 | tee $LOGDIR/data-download.log

# Step 3 - download reference and processing

bash ./ref-download.sh 2>&1 | tee $LOGDIR/ref-download.log

# Step 4 - predict putative ORFs

bash ./orf-prediction.sh 2>&1 | tee $LOGDIR/orf-prediction.log

# Step 5 - process ribo-seq data

bash ./riboseq-process.sh 2>&1 | tee $LOGDIR/riboseq-process.log

# Step 6 - process rna-seq data

bash ./rnaseq-process.sh 2>&1 | tee $LOGDIR/rnaseq-process.log

# Step 7 -  orf calling

bash ./orf-calling.sh 2>&1 | tee $LOGDIR/orf-calling.log

echo "Finished. Please find results in ./out/calling/info_table."
