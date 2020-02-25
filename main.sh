#!/bin/bash

# Master script

# Parse argument, REF - https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -d|--outdir)
    OUTDIR="$2"
    shift # past argument
    shift # past value
    ;;
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

#OUTDIRPATH=`readlink -e $OUTDIR`
#LOGDIR=$OUTDIRPATH/log

LOGDIR=$OUTDIR/log

echo $OUTDIR
echo $OUTDIRPATH
echo $LOGDIR

#mkdir -p $OUTDIRPATH
mkdir -p $OUTDIR

# Exit if received a permission denied message
if finderrors=$(! find "$LOGDIR" -depth -type d 2>&1 1>/dev/null)
    then
        echo
        echo "No permission to create the specified directory: $OUTDIRPATH"
        echo "Please choose a different location for your project."
        exit $?
fi

mkdir -p $LOGDIR
echo "Output messages can be found in $LOGDIR"
echo

exit 0

# Step 1 - check dependencies
bash ./module-check.sh 2>&1 | tee $LOGDIR/module-check.log

# Step 2 - download sample data

bash ./data-download.sh 2>&1 | tee $LOGDIR/data-download.log

# Step 3 - download reference and processing

bash ./ref-download.sh -o $ORGANISM -r $RELEASE -t $THREAD $2>&1 | tee $LOGDIR/ref-download.log

# Step 4 - predict putative ORFs

bash ./orf-prediction.sh 2>&1 | tee $LOGDIR/orf-prediction.log

# Step 5 - process ribo-seq data

bash ./riboseq-process.sh 2>&1 | tee $LOGDIR/riboseq-process.log

# Step 6 - process rna-seq data

bash ./rnaseq-process.sh 2>&1 | tee $LOGDIR/rnaseq-process.log

# Step 7 -  orf calling

bash ./orf-calling.sh 2>&1 | tee $LOGDIR/orf-calling.log

echo "Finished. Please find results in ./out/calling/info_table."
