# This script has commands to download RNA-seq and Ribo-seq data for orf-discovery project
# B - B Cell
# T - T Cell
# BMDM - Bone Marrow-Derived Macrophages

mkdir -p B/rna-seq/sasha/LPS/fastq
mkdir -p B/ribo-seq/sasha/LPS/fastq
mkdir -p B/ribo-seq/manuel/LPS/fastq
mkdir -p B/rna-seq/manuel/LPS/fastq
mkdir -p B/ribo-seq/manuel/Resting/fastq
mkdir -p B/rna-seq/manuel/Resting/fastq

mkdir -p T/ribo-seq/seb/Activated/fastq

# Download from sierra (http://www.bioinformatics.babraham.ac.uk/sierra/sierra.pl)
# THIS STEP HAS TO BE MANUAL AT THE MOMEN
# Data location refer to Google Drive Spreadsheet

# Merge reads from different lanes on flowcell
# ref - https://www.biostars.org/p/81924/

cd B/ribo-seq/sasha/LPS/fastq
parallel --dry-run 'cat B_ribo-seq_sasha_LPS_WT{}_1.fastq.gz B_ribo-seq_sasha_LPS_WT{}_2.fastq.gz >B_ribo-seq_sasha_LPS_WT{}.fastq.gz' ::: 1 2 6 4 5
rm *_1.fastq.gz
rm *_2.fastq.gz
