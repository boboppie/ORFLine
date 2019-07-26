#!/bin/bash

# Download GENCODE reference genome and annotations

# Note: GENCODE changed its server from Sanger to EBI - ftp://ftp.ebi.ac.uk/

echo
echo "--------------------------------------------------------------------------"
echo "This script will download the reference files from GENCODE and and prepare" 
echo "them for the pipeline."
echo
echo "Prerequisite programmes:"
echo "wget"
echo "gunzip"
echo "Bowtie"
echo "STAR"
echo "HTSlib"
echo
echo "Arguments:"
echo "-o|--organism : e.g. mouse or human, default: mouse"
echo "-r|--release : GENCODE annotation release, default: M22"
echo "-t|--thread : number of computer cores, default: 4"
echo
echo "Example command: ref-download.sh -o mouse -r M22 -t 4"
echo "--------------------------------------------------------------------------"
echo

echo
echo "--------------------------------------------------------------------------"
echo "Checking the availability of tools needed..."
echo "--------------------------------------------------------------------------"
echo

 
# Declare a string array with type
declare -a CmdArray=("wget" "gunzip" "bgzip" "samtools" "tabix" "bowtie-build" "STAR")
 
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
echo "--------------------------------------------------------------------------"
echo "Checking the EBI FTP server is up running..."
echo "--------------------------------------------------------------------------"
echo

FTPSTATUS=$(wget "ftp://ftp.ebi.ac.uk/" --timeout 30 -O - 2>/dev/null | grep "pub" || echo "The site is down")

if [[ $FTPSTATUS == *"down"* ]]; then
    echo
    echo "Error: can not access to the EBI FTP server...please check the network connection." >&2
    exit 1
  else
    echo "The EBI FTP server is up."
    echo
fi


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
    ORGANISM="mouse"
fi

if [ "$RELEASE" == "" ]
then
    RELEASE=M22
fi

if [ "$THREAD" == "" ]
then
    THREAD=4
fi

echo "Input:"
echo "Orgamism - $ORGANISM"
echo "Release - $RELEASE"
echo "Thread - $THREAD"
echo


CODEBASE=.
REFDIR=./out/ref
GENCODEDIR=$REFDIR/$ORGANISM/GENCODE/$RELEASE
mkdir -p $GENCODEDIR
rm -rf $GENCODEDIR/*

### Sequence fasta

# Genome PRI (CHR + scaffolds)

echo "Downloading reference genome sequence, primary assembly..."

mkdir -p $GENCODEDIR/fasta/genome/PRI
wget -nv --tries=10 -O $GENCODEDIR/fasta/genome/PRI/genome.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_$RELEASE/GRCm38.primary_assembly.genome.fa.gz # genome assembly version may change

# Transcripts CHR

echo "Downloading transcript sequences..."

mkdir -p $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts
wget -nv --tries=10 -O $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.transcripts.fa.gz

echo "Downloading protein-coding transcript sequences..."

mkdir -p $GENCODEDIR/fasta/transcriptome/CHR/protein-coding-transcripts
wget -nv --tries=10 -O $GENCODEDIR/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.pc_transcripts.fa.gz

echo "Downloading long non-coding RNA transcript sequences..."

mkdir -p $GENCODEDIR/fasta/transcriptome/CHR/lncRNA-transcripts
wget -nv --tries=10 -O $GENCODEDIR/fasta/transcriptome/CHR/lncRNA-transcripts/lncRNA_transcripts.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.lncRNA_transcripts.fa.gz

### Annitation GTF/GFF3

# Comprehensive gene annotation CHR

echo "Downloading comprehensive gene annotation..."

mkdir -p $GENCODEDIR/annotation/CHR/comprehensive
wget -nv --tries=10 -O $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.annotation.gtf.gz
wget -nv --tries=10 -O $GENCODEDIR/annotation/CHR/comprehensive/annotation.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.annotation.gff3.gz

gunzip $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf.gz

# Long non-coding RNA gene annotation CHR

#echo "Downloading long non-coding RNA gene annotation..."

#mkdir -p $GENCODEDIR/annotation/CHR/lncRNA
#wget -O $GENCODEDIR/annotation/CHR/lncRNA/lncRNA.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.long_noncoding_RNAs.gtf.gz
#wget -O $GENCODEDIR/annotation/CHR/lncRNA/lncRNA.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$ORGANISM/release_$RELEASE/gencode.v$RELEASE.long_noncoding_RNAs.gff3.gz

wget -nv --tries=10 -O $GENCODEDIR/annotation/CHR/comprehensive/MGI-ENSEMBL.rpt http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt

#---------------------------------------
# Processing reference data

# Extract chromosomes from primary assambly 

echo "Extracting chromosomes from primary assambly..."
echo 

mkdir -p $GENCODEDIR/fasta/genome/CHR/
$CODEBASE/util/fastagrep.sh 'chr' <(gzip -dc $GENCODEDIR/fasta/genome/PRI/genome.fa.gz) >$GENCODEDIR/fasta/genome/CHR/genome.fa

# Filter for transcriptome

echo "Filtering for transcriptome..."
echo "The following biotypes will be filtered out:"
echo "IG_*, TR_*, miRNA, misc_RNA, Mt_*, rRNA, scaRNA, scRNA, snoRNA, snRNA, sRNA, ribozyme, nonsense_mediated_decay, non_stop_decay"
echo 

$CODEBASE/util/fastagrep.sh -v 'IG_*|TR_*|miRNA|misc_RNA|Mt_*|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|ribozyme|nonsense_mediated_decay|non_stop_decay' <(gzip -dc $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz) >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered.fa

grep ">" $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered.fa | sed 's/>//g' | cut -f1,7 -d"|" | tr "|" "\t" | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_Id_len.tsv

grep ">" $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered.fa | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_header.txt

grep protein_coding $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_header.txt | cut -f1 -d"|" | cut -f2 -d">" | sort | uniq >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_Ids.txt

grep ">" <(gzip -dc $GENCODEDIR/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.gz) >$GENCODEDIR/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts_header.txt

LC_ALL=C fgrep -f $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_Ids.txt $GENCODEDIR/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts_header.txt | cut -c 2- | tr "|" "\t" | awk -v OFS='\t' '{print $1, $8, $9, $10}' >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo.tsv

awk '$2~/CDS/ && $3==""' $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo.tsv | awk -v OFS='\t' '{print $1, "UTR5:NULL-NULL", $2, "UTR3:NULL-NULL"}' >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_CDSonly.tsv 
awk '$2~/CDS/ && $3~/UTR3/' $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo.tsv | awk -v OFS='\t' '{print $1, "UTR5:NULL-NULL", $2, $3}' >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_CDSandUTR3.tsv
awk '$3~/CDS/ && $2~/UTR5/ && $4==""' $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo.tsv | awk -v OFS='\t' '{print $1, $2, $3, "UTR3:NULL-NULL"}' >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_UTR5andCDS.tsv
awk '$3~/CDS/ && $2~/UTR5/ && $4~/UTR3/' $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo.tsv >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_UTR5andCDSandUTR3.tsv
cat $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_*.tsv | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted.tsv
sed "s/UTR5://" $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted.tsv | sed "s/CDS://" | sed "s/UTR3://" | tr "-" "\t" | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted_simplified.tsv

# Stop codon inclueded in CDS
cut -f1,3 $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted.tsv | tr "-" "\t" | sed "s/CDS://" | awk -v OFS=":" '{if ($3-$2+1 <= 303) print $1,$2,$3}' | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_smORFs_canonical.txt
cut -f1,3 $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted.tsv | tr "-" "\t" | sed "s/CDS://" | awk -v OFS=":" '{if ($3-$2+1 > 303) print $1,$2,$3}' | sort >$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_non-smORFs_canonical.txt


echo "Creating filtered protein-coding annotation..."
echo

cat $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf | awk '{if($18=="\"protein_coding\";" && $0~"level (1|2);" && $0!~"tag \"seleno\";" && $0!~"cds_end_NF" && $0!~"cds_start_NF"){print $0}}' | sort -k1,1 -k4,4n | bgzip >$GENCODEDIR/annotation/CHR/comprehensive/pc_level1_2_noSeleno.gtf.gz
tabix -p gff $GENCODEDIR/annotation/CHR/comprehensive/pc_level1_2_noSeleno.gtf.gz

awk '{if ($3 == "CDS" || $3 == "stop_codon") print}' $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf >$GENCODEDIR/annotation/CHR/comprehensive/annotation_CDS_with_STOP.gtf
awk '{if($18=="\"protein_coding\";"){print $0}}' $GENCODEDIR/annotation/CHR/comprehensive/annotation_CDS_with_STOP.gtf | sort -k1,1 -k4,4n | bgzip >$GENCODEDIR/annotation/CHR/comprehensive/CDS_PC_sorted.gtf.gz
tabix -p gff $GENCODEDIR/annotation/CHR/comprehensive/CDS_PC_sorted.gtf.gz

awk '$3=="transcript"' $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf | awk '{print $12}' | cut -f2 -d'"' | sort >$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids.txt

cut -f1 -d. $GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids.txt >$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids_noVersion.txt
paste $GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids_noVersion.txt $GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids.txt | sort >$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids.tsv

# Extract gene description from MGI
# Name the result to All_transcript_ids_geneName_description.tsv
grep -v '^#' $GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf | cut -f9 | grep transcript_id | sed "s/\"//g; s/;//g" | cut -f2,4,8 -d " " | sed "s/ /\t/g" | cut -f1,3,5 | sort | uniq >$GENCODEDIR/annotation/CHR/comprehensive/geneId_ver-txId_ver-symbol.tsv
paste <(cut -f1 $GENCODEDIR/annotation/CHR/comprehensive/geneId_ver-txId_ver-symbol.tsv | sed 's/\./\t/g') <(cut -f2-3 $GENCODEDIR/annotation/CHR/comprehensive/geneId_ver-txId_ver-symbol.tsv) >$GENCODEDIR/annotation/CHR/comprehensive/geneId-txId_ver-symbol.tsv
cut -f2,3,6 $GENCODEDIR/annotation/CHR/comprehensive/MGI-ENSEMBL.rpt | awk -v OFS='\t' -F '\t' '{print $3,$1,$2}' | sort | uniq >$GENCODEDIR/annotation/CHR/comprehensive/MGI-ENSEMBL-geneId-symbol-description.tsv
join -a1 -a2 -t $'\t' $GENCODEDIR/annotation/CHR/comprehensive/geneId-txId_ver-symbol.tsv $GENCODEDIR/annotation/CHR/comprehensive/MGI-ENSEMBL-geneId-symbol-description.tsv | cut -f2,3,5 | sort | uniq >$GENCODEDIR/annotation/CHR/comprehensive/txId_ver-symbol-description.tsv

sort $CODEBASE/ref/annotation/All_transcript_ids_geneName_description.tsv >$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids_geneName_description_sorted.tsv
join -t $'\t' -a1 -a2 $GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids.tsv $GENCODEDIR/annotation/CHR/comprehensive/All_transcript_ids_geneName_description_sorted.tsv | cut -f2-4 >$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_idsWithVersion_geneName_description.tsv

# Contaimination sequences

# Some sequences have pre-downloaded to the code repository
echo "Merging rRNA/tRNA sequences..."
echo

mkdir -p $GENCODEDIR/contaimination/bowtie_index

$CODEBASE/util/fastagrep.sh 'rRNA|Mt_rRNA|Mt_tRNA|snRNA|snoRNA|misc_RNA|miRNA' <(gzip -dc $GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz) > $GENCODEDIR/contaimination/GENCODE_contaimination.fa

cat $GENCODEDIR/contaimination/GENCODE_contaimination.fa $CODEBASE/ref/sequences/UCSC_${ORGANISM}_*.fa $CODEBASE/ref/sequences/NCBI_${ORGANISM}_rRNA.fa >$GENCODEDIR/contaimination/contaimination.fa


# softlinking to transcript sequences and annotation

echo "Making reference files global (in ref/pipeline_global)..."
echo
mkdir -p $REFDIR/pipeline_global

# remove old soft links
if [ "$(ls -A $REFDIR/pipeline_global)" ]; then
     rm -rf $REFDIR/pipeline_global/*
fi

ln -s ../../../$GENCODEDIR/fasta/genome/CHR/genome.fa $REFDIR/pipeline_global/genome.fa
ln -s ../../../$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered.fa $REFDIR/pipeline_global/transcriptome.fa
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/annotation.gtf $REFDIR/pipeline_global/annotation.gtf
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/annotation.gff3.gz $REFDIR/pipeline_global/annotation.gff3.gz
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/pc_level1_2_noSeleno.gtf.gz $REFDIR/pipeline_global/pc.gtf.gz
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/pc_level1_2_noSeleno.gtf.gz.tbi $REFDIR/pipeline_global/pc.gtf.gz.tbi
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/annotation_CDS_with_STOP.gtf $REFDIR/pipeline_global/annotation_CDS_with_STOP.gtf
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/CDS_PC_sorted.gtf.gz $REFDIR/pipeline_global/CDS_PC_sorted.gtf.gz
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/CDS_PC_sorted.gtf.gz.tbi $REFDIR/pipeline_global/CDS_PC_sorted.gtf.gz.tbi
ln -s ../../../$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_smORFs_canonical.txt $REFDIR/pipeline_global/transcripts_biotype_filtered_pc_posInfo_smORFs_canonical.txt
ln -s ../../../$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_pc_posInfo_formatted_simplified.tsv $REFDIR/pipeline_global/transcripts_biotype_filtered_pc_posInfo_formatted_simplified.tsv
ln -s ../../../$GENCODEDIR/fasta/transcriptome/CHR/all-transcripts/transcripts_biotype_filtered_header.txt $REFDIR/pipeline_global/transcripts_biotype_filtered_header.txt
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/All_transcript_idsWithVersion_geneName_description.tsv $REFDIR/pipeline_global/All_transcript_idsWithVersion_geneName_description.tsv
ln -s ../../../$GENCODEDIR/annotation/CHR/comprehensive/txId_ver-symbol-description.tsv $REFDIR/pipeline_global/txId_ver-symbol-description.tsv

# Building index
echo "Building Bowtie index..."
echo

bowtie-build $GENCODEDIR/contaimination/contaimination.fa $GENCODEDIR/contaimination/bowtie_index

ln -s ../../../$GENCODEDIR/contaimination $REFDIR/pipeline_global/contaimination

# STAR indexing

echo "Generateing STAR index for 50/100bp read..."
echo

SJ_GTF=$REFDIR/pipeline_global/annotation.gtf

mkdir -p $GENCODEDIR/star/50bp
STAR --runThreadN $THREAD --runMode genomeGenerate --genomeDir $GENCODEDIR/star/50bp --genomeFastaFiles $GENCODEDIR/fasta/genome/CHR/genome.fa --sjdbGTFfile $SJ_GTF --sjdbOverhang 49 

#mkdir -p $REFDIR/$ORGANISM/star/75bp
#STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $REFDIR/$ORGANISM/star/75bp --genomeFastaFiles $GENCODEDIR/fasta/genome/CHR/genome.fa --sjdbGTFfile $SJ_GTF --sjdbOverhang 74 

mkdir -p $GENCODEDIR/star/100bp
STAR --runThreadN $THREAD --runMode genomeGenerate --genomeDir $GENCODEDIR/star/100bp --genomeFastaFiles $GENCODEDIR/fasta/genome/CHR/genome.fa --sjdbGTFfile $SJ_GTF --sjdbOverhang 99

ln -s ../../../$GENCODEDIR/star $REFDIR/pipeline_global/star

