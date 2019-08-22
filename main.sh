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

# Step 1 - check dependencies
bash ./module-check.sh >$LOGDIR/module-check.log 2>&1

# Step 2 - download sample data

bash ./data-download.sh >$LOGDIR/data-download.log 2>&1

# Step 3 - download reference and processing

bash ./ref-download.sh >$LOGDIR/ref-download.log 2>&1

# Step 4 - predict putative ORFs

bash ./orf-prediction.sh >$LOGDIR/orf-prediction.log 2>&1

# Step 5 - process ribo-seq data

bash ./riboseq-process.sh >$LOGDIR/riboseq-process.log 2>&1

# Step 6 - process rna-seq data

bash ./rnaseq-process.sh >$LOGDIR/rnaseq-process.log 2>&1

# Step 7 -  orf calling

bash ./orf-calling.sh >$LOGDIR/orf-calling.log 2>&1

echo "Finished. Please find results in ./out/calling/info_table."
