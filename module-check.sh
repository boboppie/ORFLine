#!/bin/bash

echo
echo "The following tools or packages are needed to run the pipeline:"
echo
echo "----------------------------------------------------------------------------------"
echo "Samtools and HTSlib"
echo "bedtools"
echo "BEDOPS"
echo "Bowtie"
echo "STAR"
echo "FastQC"
echo "Trim Galore"
echo "plastid"
echo "StringTie"
echo "sra tools"
echo "EMBOSS"
echo "GNU Parallel"
echo "R"
echo
echo "R/Bioconductor pakcages:"
echo
echo "riboSeqR"
echo "GenomicFeatures"
echo "rtracklayer"
echo "----------------------------------------------------------------------------------"
echo


echo
echo "----------------------------------------------------------------------------------"
echo "Checking the availability of tools needed for the pipeline..."
echo "----------------------------------------------------------------------------------"
echo

 
# Declare a string array with type
declare -a CmdArray=("samtools" "tabix" "bedtools" "bedops" "bowtie" "fastqc" "fastq-dump" "trim_galore" "psite" "STAR" "stringtie" "transeq" "parallel" "R")
 
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

echo
echo "----------------------------------------------------------------------------------"
echo "Checking the availability of Bioconductor packages needed for the pipeline..."
echo "----------------------------------------------------------------------------------"
echo

declare -a RpkgArray=("rtracklayer" "GenomicFeatures" "riboSeqR")

for pkg in "${RpkgArray[@]}"; do
  RESULT=$(Rscript <(echo "\"$pkg\" %in% rownames(installed.packages())"))

  if [[ $RESULT == *"FALSE"* ]]; then
    echo
    echo "Error: Bioconductor package $pkg is not installed." >&2
    exit 1
  else
    echo "Bioconductor package $pkg is available."  
  fi
done

echo "All installed."
