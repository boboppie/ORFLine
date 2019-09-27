#!/bin/bash

echo
echo "--------------------------------------------------------------------------------------------"
echo "This script will call actively translated smORFs using Ribo-Seq and RNA-Seq BAMs."
echo
echo "Arguments:"
echo "-o|--organism : e.g. mouse or human, default: mouse"
echo "-x|--taxid : taxid of the organism, e.g. 10090 for mouse, default: 10090"
echo "-m|--maxl : max length of ribosome protected fragment, default: 32"
echo "-n|--minl: min length of ribosome protected fragment, default: 28"
echo "-t|--thread : number of computer cores, default: 4"
echo
echo "Example command: orf-calling.sh -o mouse -x 10090 -m 32 -n 28 -t 8"
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
    -x|--taxid)
    TAXID="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--maxl)
    MAXL="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--minl)
    MINL="$2"
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

if [ "$TAXID" == "" ]
then
    TAXID=10090
fi

if [ "$MAXL" == "" ]
then
    MAXL=32
fi

if [ "$MINL" == "" ]
then
    MINL=28
fi

if [ "$THREAD" == "" ]
then
    THREAD=4
fi


echo "Input:"
echo "Orgamism - $ORGANISM"
echo "Taxid - $TAXID"
echo "Max length - $MAXL"
echo "Min length - $MINL"
echo "Thread - $THREAD"
echo

OUTPATH=./out/calling
mkdir -p $OUTPATH

CODEBASE=.

RIBOBAMDIR=./out/process/ribo/star-genome
RNABAMDIR=./out/process/rna/star-genome
RNA_EXPRESSION_DIR=./out/process/rna/stringtie

PSITEOFFSET=./out/process/ribo/plastid/psite/*_star_genome_p_offsets.txt
GTF=./out/ref/pipeline_global/annotation.gtf
GFF=./out/ref/pipeline_global/annotation.gff3.gz
CDSGTF=./out/ref/pipeline_global/annotation_CDS_with_STOP.gtf
CDS_PC_GTF=./out/ref/pipeline_global/CDS_PC_sorted.gtf.gz
TRANSCRIPTOME=./out/ref/pipeline_global/transcriptome.fa
REFGENOME=./out/ref/pipeline_global/genome.fa
GENEDESCRIPTION=./out/ref/pipeline_global/txId_ver-symbol-description.tsv


# Step 1:
# BAM filter: only keep RPF in given range [MINL, MAXL]

echo "Start BAM filter..."

mkdir -p $OUTPATH/bam_filter

# Samtools to filter reads
samtools view -h $RIBOBAMDIR/*_star_genome_q255.sorted.bam | awk -v minl="$MINL" -v maxl="$MAXL" 'length($10) >= minl && length($10) <= maxl || $1 ~ /^@/' | samtools view -bS - > $OUTPATH/bam_filter/filtered_${MINL}-${MAXL}nt.bam
samtools index $OUTPATH/bam_filter/filtered_${MINL}-${MAXL}nt.bam
BAM=$OUTPATH/bam_filter/filtered_${MINL}-${MAXL}nt.bam

echo "Finished BAM filter."
echo

# Step 2: 
# Filter out smORFs having 0 read mapped to the region

echo "Start smORF read filter..."

mkdir -p $OUTPATH/counts_bedtools

# Need to use -split as the input is BED12
parallel -j $THREAD "bedtools coverage -counts -split -a ./out/prediction/orfs_{}.smORFs.bed -b $BAM >$OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.counts" ::: ATG CTG TTG GTG

# Only filter out 0 read ORFs (can not set a cutoff by looking at the distribution)
parallel -j $THREAD "awk '\$13 > 0' $OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.counts >$OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.greaterThanZeroReads.counts" ::: ATG CTG TTG GTG

parallel -j $THREAD "cut -f1-12 $OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.greaterThanZeroReads.counts >$OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.greaterThanZeroReads.bed" ::: ATG CTG TTG GTG

echo "Finished smORF read filter."
echo

# Step 3: 
# Filter out smORFs having 0 RPF mapped to the region

echo "Start smORF RPF filter..."

mkdir -p $OUTPATH/counts_plastid

parallel -j $THREAD "python $CODEBASE/util/get_count_vectors_updated.py --annotation_files $OUTPATH/counts_bedtools/orfs_{}.smORFs.bedtools.greaterThanZeroReads.bed --annotation_format BED --count_files $BAM --min_length $MINL --max_length $MAXL --fiveprime_variable --offset $PSITEOFFSET --out_prefix {} $OUTPATH/counts_plastid" ::: ATG CTG TTG GTG


parallel -j $THREAD "awk '{for(i=2;i<=NF;i++) t+=\$i; print \$1\"\t\"t; t=0}' $OUTPATH/counts_plastid/ORFs_{}.txt >$OUTPATH/counts_plastid/ORFs_{}_rowSum.txt" ::: ATG CTG TTG GTG

cat $OUTPATH/counts_plastid/ORFs_*TG_rowSum.txt >$OUTPATH/counts_plastid/ORFs_ALL_rowSum.txt

# Set cutoff to 0
parallel -j $THREAD "paste $OUTPATH/counts_plastid/ORFs_{}_rowSum.txt $OUTPATH/counts_plastid/ORFs_{}.txt | awk '{if (\$2 > 0) print}' | cut -f3- >$OUTPATH/counts_plastid/ORFs_{}_translated.txt; paste $OUTPATH/counts_plastid/ORFs_{}_rowSum.txt $OUTPATH/counts_plastid/ORFs_{}.txt | awk '{if (\$2 == 0) print}' | cut -f3- >$OUTPATH/counts_plastid/ORFs_{}_nonTranslated.txt" ::: ATG CTG TTG GTG

echo "Finished smORF RPF filter."
echo


# Step 4: 
# Filter out lowly or non expressed host transcripts

echo "Start transcript expression filter..."

mkdir -p $OUTPATH/tx_expressed_filter

parallel -j $THREAD "LC_ALL=C fgrep -f ./out/process/rna/stringtie/active_transcript_Ids.txt $OUTPATH/counts_plastid/ORFs_{}_translated.txt >$OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed.txt" ::: ATG CTG TTG GTG

parallel -j $THREAD "cut -f1 $OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed.txt >$OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed_ORFId.txt" ::: ATG CTG TTG GTG

parallel -j $THREAD "LC_ALL=C fgrep -f $OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed_ORFId.txt ./out/prediction/orfs_{}.smORFs.bed >$OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed.bed" ::: ATG CTG TTG GTG

cat $OUTPATH/tx_expressed_filter/ORFs_*_translated_expressed.bed >$OUTPATH/tx_expressed_filter/ORFs_ALL_translated_expressed.bed

echo "Finished transcript expression filter."
echo


# Step 5: 
# Assign labels to smORFs

echo "Start label filter..."

mkdir -p $OUTPATH/label_filter

awk -v OFS='\t' '{print $1"_"$2"_"$3"_"$6"_"$10"_"$11"_"$12, $4}' $OUTPATH/tx_expressed_filter/ORFs_ALL_translated_expressed.bed >$OUTPATH/label_filter/regions.tsv

cut -f2 $OUTPATH/label_filter/regions.tsv | sort >$OUTPATH/label_filter/regions_ORFIds.txt

awk -v FS=':' -v OFS="\t" '{print $1":"$4":"$5, $1":"$2":"$3":"$4":"$5":"$6":"$7}' $OUTPATH/label_filter/regions_ORFIds.txt | sort >$OUTPATH/label_filter/regions_ORFIds_extra.tsv

LC_ALL=C fgrep -f ./out/ref/pipeline_global/transcripts_biotype_filtered_pc_posInfo_smORFs_canonical.txt $OUTPATH/label_filter/regions_ORFIds_extra.tsv | cut -f2 | sort >$OUTPATH/label_filter/regions_ORFIds_canonical.txt

LC_ALL=C fgrep -v -f $OUTPATH/label_filter/regions_ORFIds_canonical.txt $OUTPATH/label_filter/regions_ORFIds.txt | sort >$OUTPATH/label_filter/regions_ORFIds_nonCanonical.txt

grep protein_coding $OUTPATH/label_filter/regions_ORFIds_nonCanonical.txt >$OUTPATH/label_filter/regions_ORFIds_nonCanonical_pc.txt 
grep -v protein_coding $OUTPATH/label_filter/regions_ORFIds_nonCanonical.txt >$OUTPATH/label_filter/regions_ORFIds_nonCanonical_noncoding.txt 

awk '{print $1, "canonical"}' $OUTPATH/label_filter/regions_ORFIds_canonical.txt | cut -f2 -d' ' >$OUTPATH/label_filter/regions_ORFIds_canonical_label.txt
awk '{print $1, "unsigned"}' $OUTPATH/label_filter/regions_ORFIds_nonCanonical_pc.txt | cut -f2 -d' ' >$OUTPATH/label_filter/regions_ORFIds_nonCanonical_pc_label.txt
awk '{print $1, "noncoding"}' $OUTPATH/label_filter/regions_ORFIds_nonCanonical_noncoding.txt | cut -f2 -d' ' >$OUTPATH/label_filter/regions_ORFIds_nonCanonical_noncoding_label.txt 

# Sign label to "unsigned protein coding ORFs"
cut -f1,4,5 -d: $OUTPATH/label_filter/regions_ORFIds_nonCanonical_pc.txt | tr ":" "\t" | paste - $OUTPATH/label_filter/regions_ORFIds_nonCanonical_pc.txt | join ./out/ref/pipeline_global/transcripts_biotype_filtered_pc_posInfo_formatted_simplified.tsv - | tr " " "\t" >$OUTPATH/label_filter/label_calling_info.tsv

awk -v OFS="\t" '{if($8 >= $2 && $9 <= $3) print $10, "five_prime"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_five_prime.tsv
awk -v OFS="\t" '{if($8 >= $6 && $9 <= $7) print $10, "three_prime"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_three_prime.tsv
awk -v OFS="\t" '{if($8 > $4 && $9 < $5) print $10, "within"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_within.tsv
awk -v OFS="\t" '{if($8 > $4 && $9 == $5) print $10, "canonical_truncated"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_canonical_truncated.tsv
awk -v OFS="\t" '{if($8 < $4 && $9 == $5) print $10, "canonical_extented"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_canonical_extented.tsv
awk -v OFS="\t" '{if($8 >= $2 && $8 <= $3 && $9 >= $4 && $9 < $5) print $10, "five_prime_overlap"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_five_prime_overlap.tsv
awk -v OFS="\t" '{if($8 > $4 && $8 <= $5 && $9 >= $6 && $9 <= $7) print $10, "three_prime_overlap"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_three_prime_overlap.tsv
awk -v OFS="\t" '{if($8 == $4 && $9 < $5) print $10, "seleno"}' $OUTPATH/label_filter/label_calling_info.tsv >$OUTPATH/label_filter/label_calling_out_seleno.tsv

cat $OUTPATH/label_filter/label_calling_out_* | cut -f1 | sort >$OUTPATH/label_filter/label_calling_out_ORFId.txt
LC_ALL=C fgrep -v -f $OUTPATH/label_filter/label_calling_out_ORFId.txt $OUTPATH/label_filter/label_calling_info.tsv | cut -f10 | awk -v OFS="\t" '{print $1, "uncertain"}' >$OUTPATH/label_filter/label_calling_out_uncertain.tsv
cat $OUTPATH/label_filter/label_calling_out_five_prime.tsv $OUTPATH/label_filter/label_calling_out_three_prime.tsv $OUTPATH/label_filter/label_calling_out_within.tsv $OUTPATH/label_filter/label_calling_out_canonical_truncated.tsv $OUTPATH/label_filter/label_calling_out_canonical_extented.tsv $OUTPATH/label_filter/label_calling_out_five_prime_overlap.tsv $OUTPATH/label_filter/label_calling_out_three_prime_overlap.tsv $OUTPATH/label_filter/label_calling_out_uncertain.tsv >$OUTPATH/label_filter/label_calling_out_final.tsv
sort $OUTPATH/label_filter/label_calling_out_final.tsv >$OUTPATH/label_filter/label_calling_out_final_sorted.tsv
cut -f1 $OUTPATH/label_filter/label_calling_out_final_sorted.tsv >$OUTPATH/label_filter/unsigned_ORFId.txt
cut -f2 $OUTPATH/label_filter/label_calling_out_final_sorted.tsv >$OUTPATH/label_filter/unsigned_label.txt

LC_ALL=C fgrep -f $OUTPATH/label_filter/regions_ORFIds_canonical.txt $OUTPATH/tx_expressed_filter/ORFs_ALL_translated_expressed.bed >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_canonical.bed
LC_ALL=C fgrep -f $OUTPATH/label_filter/unsigned_ORFId.txt $OUTPATH/tx_expressed_filter/ORFs_ALL_translated_expressed.bed >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_nonCanonical_pc.bed
LC_ALL=C fgrep -f $OUTPATH/label_filter/regions_ORFIds_nonCanonical_noncoding.txt $OUTPATH/tx_expressed_filter/ORFs_ALL_translated_expressed.bed >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_nonCanonical_noncoding.bed

paste $OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_canonical.bed $OUTPATH/label_filter/regions_ORFIds_canonical_label.txt >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12Plus_canonical_label.bed
paste <(sort -k4 $OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_nonCanonical_pc.bed) $OUTPATH/label_filter/unsigned_label.txt >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12Plus_nonCanonical_pc_label.bed
paste $OUTPATH/label_filter/ORFScore_ALL_filtered_BED12_nonCanonical_noncoding.bed $OUTPATH/label_filter/regions_ORFIds_nonCanonical_noncoding_label.txt >$OUTPATH/label_filter/ORFScore_ALL_filtered_BED12Plus_nonCanonical_noncoding_label.bed

cat $OUTPATH/label_filter/ORFScore_ALL_filtered_BED12Plus_*_label.bed >$OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed

cut -f4 $OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed >$OUTPATH/label_filter/ORFId_label_filter.txt

echo "Finished label filter."
echo


# Step 6:
# Calulate ORFScore

echo "Start ORFScore filter..."

mkdir -p $OUTPATH/orfscore_filter

parallel -j $THREAD "grep {} $OUTPATH/label_filter/ORFId_label_filter.txt >$OUTPATH/orfscore_filter/ORFId_label_filter_{}.txt" ::: ATG CTG TTG GTG
parallel -j $THREAD "LC_ALL=C fgrep -f $OUTPATH/orfscore_filter/ORFId_label_filter_{}.txt $OUTPATH/tx_expressed_filter/ORFs_{}_translated_expressed.txt >$OUTPATH/orfscore_filter/ORFs_{}_translated_expressed.txt" ::: ATG CTG TTG GTG

parallel -j $THREAD "Rscript $CODEBASE/util/ORFScore_updated.R $OUTPATH/orfscore_filter/ORFs_{}_translated_expressed.txt $OUTPATH/orfscore_filter {}" ::: ATG CTG TTG GTG

cat $OUTPATH/orfscore_filter/ORFScore_*TG.tsv | sort -k2 -g | awk '{if($2!="NA") print}' >$OUTPATH/orfscore_filter/ORFScore_all_noNA.tsv
cat $OUTPATH/orfscore_filter/ORFScore_all_noNA.tsv | awk '{if($2 >0 && $6 >= 0.1) print}' >$OUTPATH/orfscore_filter/ORFScore_all_filteredByCov.tsv

cut -f1 $OUTPATH/orfscore_filter/ORFScore_all_filteredByCov.tsv >$OUTPATH/orfscore_filter/ORFId_orfscore_filter.txt

echo "Finished ORFScore filter."
echo


# Step 7:
# Region filter: filter out prediction overlapping with exons

echo "Start region filter..."

mkdir -p $OUTPATH/region_filter

LC_ALL=C fgrep -f $OUTPATH/orfscore_filter/ORFId_orfscore_filter.txt $OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed >$OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed
cut -f1-12 $OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed >$OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED12.bed
bed12ToBed6 -i $OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED12.bed >$OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED6.bed

Rscript $CODEBASE/util/region_filter.R $OUTPATH/region_filter ../../../$CDSGTF ../../../$BAM $THREAD

cut -f6 $OUTPATH/region_filter/orfs_readRatio_canonical.tsv | sort >$OUTPATH/region_filter/filtered_canonical_ORFId.txt

# Kepp all $15 == "NA" or $15 < 1
LC_ALL=C fgrep -v -f $OUTPATH/region_filter/filtered_canonical_ORFId.txt $OUTPATH/region_filter/orfs_readRatio.tsv | awk -F $'\t' '$15 == "NA" || $15 < 1' | cut -f6 >$OUTPATH/region_filter/filtered_NA_ORFId.txt
cat $OUTPATH/region_filter/filtered_canonical_ORFId.txt $OUTPATH/region_filter/filtered_NA_ORFId.txt | sort >$OUTPATH/region_filter/ORFId_region_filter.txt

echo "Finished region filter."
echo


# Step 8:
# Calculate RRS 

echo "Start rrs filter..."

mkdir -p $OUTPATH/rrs_filter

LC_ALL=C fgrep -f $OUTPATH/region_filter/ORFId_region_filter.txt $OUTPATH/region_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed >$OUTPATH/rrs_filter/filtered_BED12Plus.bed
bed12ToBed6 -i <(cut -f1-12 $OUTPATH/rrs_filter/filtered_BED12Plus.bed) >$OUTPATH/rrs_filter/filtered_BED6.bed

# smORF length
paste -d'\t' <(cut -f4 $OUTPATH/rrs_filter/filtered_BED12Plus.bed) <(cut -f11 $OUTPATH/rrs_filter/filtered_BED12Plus.bed | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1') >$OUTPATH/rrs_filter/filtered_smORFs_length.tsv

# Total mapped reads for Ribo and RNA-seq samples
samtools view -F 4 -c $OUTPATH/bam_filter/filtered_*.bam >$OUTPATH/rrs_filter/rp_totalMappedReads_merged.tsv

cut -f1,3 $RNA_EXPRESSION_DIR/expressed.tsv >$OUTPATH/rrs_filter/rna_fpkm_merged.tsv

grep canonical $OUTPATH/rrs_filter/filtered_BED12Plus.bed | cut -f4 >$OUTPATH/rrs_filter/ORFId_canonical.txt
cut -f1 -d: $OUTPATH/rrs_filter/ORFId_canonical.txt | sort | uniq >$OUTPATH/rrs_filter/txId_canonical.txt
LC_ALL=C fgrep -f $OUTPATH/rrs_filter/txId_canonical.txt <(gunzip -c $GFF) | awk '$3 == "three_prime_UTR"' >$OUTPATH/rrs_filter/txId_canonical_UTR3.gff3
gff2bed <$OUTPATH/rrs_filter/txId_canonical_UTR3.gff3 | cut -f1-6 | sed 's/UTR3://g' | sort -k4 >$OUTPATH/rrs_filter/txId_canonical_UTR3_BED6.bed
cut -f1 -d: $OUTPATH/rrs_filter/ORFId_canonical.txt >$OUTPATH/rrs_filter/ORFId_canonical_txId.txt
paste $OUTPATH/rrs_filter/ORFId_canonical_txId.txt $OUTPATH/rrs_filter/ORFId_canonical.txt | sort -k1 >$OUTPATH/rrs_filter/ORFId_canonical_txId.tsv
join -1 4 -2 1 $OUTPATH/rrs_filter/txId_canonical_UTR3_BED6.bed $OUTPATH/rrs_filter/ORFId_canonical_txId.tsv | awk -v OFS='\t' '{print $2, $3, $4, $7":pUTR3", $5, $6}' >$OUTPATH/rrs_filter/UTR3_canonical_BED6.bed

grep -v canonical $OUTPATH/rrs_filter/filtered_BED12Plus.bed >$OUTPATH/rrs_filter/filtered_nonCanonical.bed 
cat $OUTPATH/rrs_filter/filtered_nonCanonical.bed | cut -f4 | cut -f1 -d: | sort | uniq >$OUTPATH/rrs_filter/txId_nonCanonical.txt
LC_ALL=C fgrep -f $OUTPATH/rrs_filter/txId_nonCanonical.txt ./out/ref/pipeline_global/transcripts_biotype_filtered_header.txt >$OUTPATH/rrs_filter/txHeader_nonCanonical.txt
cut -c 2- $OUTPATH/rrs_filter/txHeader_nonCanonical.txt | xargs -n 1 samtools faidx $TRANSCRIPTOME >$OUTPATH/rrs_filter/tx_nonCanonical.fa


# Count RPF for all protein coding transcript CDSs using merged BAM
python $CODEBASE/util/get_count_vectors_updated.py --annotation_files $CDS_PC_GTF --annotation_format GTF2 --count_files $OUTPATH/bam_filter/filtered*.bam --min_length $MINL --max_length $MAXL --fiveprime_variable --offset $PSITEOFFSET --out_prefix five_prime_mainORFCDS_merged $OUTPATH/rrs_filter; awk -v OFS="\t" '{for(i=2;i<=NF;i++) t+=$i; print $1,t; t=0}' $OUTPATH/rrs_filter/ORFs_five_prime_mainORFCDS_merged.txt | sort >$OUTPATH/rrs_filter/CDS_PC_RPF_merged.tsv

gunzip -c $CDS_PC_GTF | gtf2bed | awk -v OFS="\t" '{print $1,$2,$3,$13,$5,$6}' | sed 's/"//g; s/;//g;' >$OUTPATH/rrs_filter/CDS_PC_BED6.bed 
awk -v OFS='\t' '{print $4, $3-$2}' $OUTPATH/rrs_filter/CDS_PC_BED6.bed | awk -v OFS='\t' '{a[$1] += $2} END{for (i in a) print i, a[i]}' | sort >$OUTPATH/rrs_filter/CDS_PC_length.tsv

Rscript $CODEBASE/util/preRRS.R $OUTPATH/rrs_filter ../../../$GTF $TAXID ../../../$BAM tx_nonCanonical.fa $THREAD

# for the merged
LC_ALL=C fgrep -f $OUTPATH/region_filter/ORFId_region_filter.txt $OUTPATH/counts_plastid/ORFs_ALL_rowSum.txt | sort >$OUTPATH/rrs_filter/Ribo_ALL_ORF_sum.tsv

python $CODEBASE/util/get_count_vectors_updated.py --annotation_files $OUTPATH/rrs_filter/UTR3_canonical_BED6.bed --annotation_format BED --count_files $BAM --min_length $MINL --max_length $MAXL --fiveprime_variable --offset $PSITEOFFSET --out_prefix UTR3_canonical $OUTPATH/rrs_filter
awk -v OFS="\t" '{for(i=2;i<=NF;i++) t+=$i; print $1,t; t=0}' $OUTPATH/rrs_filter/ORFs_UTR3_canonical.txt >$OUTPATH/rrs_filter/ORFs_UTR3_canonical_rowSum.tsv
awk -v OFS="\t" '{a[$1]+=$2}END{for(i in a) print i,a[i]}' $OUTPATH/rrs_filter/ORFs_UTR3_canonical_rowSum.tsv >$OUTPATH/rrs_filter/Ribo_canonical_UTR3_sum.tsv

python $CODEBASE/util/get_count_vectors_updated.py --annotation_files $OUTPATH/rrs_filter/pUTR3_BED12.bed --annotation_format BED --count_files $BAM --min_length $MINL --max_length $MAXL --fiveprime_variable --offset $PSITEOFFSET --out_prefix UTR3_nonCanonical $OUTPATH/rrs_filter
awk -v OFS="\t" '{for(i=2;i<=NF;i++) t+=$i; print $1,t; t=0}' $OUTPATH/rrs_filter/ORFs_UTR3_nonCanonical.txt >$OUTPATH/rrs_filter/Ribo_nonCanonical_UTR3_sum.tsv

cat $OUTPATH/rrs_filter/Ribo_canonical_UTR3_sum.tsv $OUTPATH/rrs_filter/Ribo_nonCanonical_UTR3_sum.tsv | sort >$OUTPATH/rrs_filter/Ribo_ALL_UTR3_sum.tsv

cut -f1-12 $OUTPATH/rrs_filter/filtered_BED12Plus.bed >$OUTPATH/rrs_filter/filtered_BED12.bed

# UTR3 length
awk -v OFS="\t" '{print $4, $3-$2}' $OUTPATH/rrs_filter/UTR3_canonical_BED6.bed | awk -v OFS="\t" '{a[$1]+=$2}END{for(i in a) print i,a[i]}' >$OUTPATH/rrs_filter/UTR3_canonical_length.tsv
awk -v OFS="\t" '{print $4, $3-$2}' $OUTPATH/rrs_filter/pUTR3_BED6.bed | awk -v OFS="\t" '{a[$1]+=$2}END{for(i in a) print i,a[i]}' >$OUTPATH/rrs_filter/pUTR3_length.tsv
cat $OUTPATH/rrs_filter/UTR3_canonical_length.tsv $OUTPATH/rrs_filter/pUTR3_length.tsv | sort >$OUTPATH/rrs_filter/UTR3_length.tsv

Rscript $CODEBASE/util/RRS_merged.R $OUTPATH/rrs_filter ../../../$RNABAMDIR/

# Tidy up translation efficiency results, remove (Inf|NA|NaN)
grep -v -E 'NA|Inf|NaN' $OUTPATH/rrs_filter/TE_smORF_teMerged.tsv >$OUTPATH/rrs_filter/TE_smORF_teMerged_tidy.tsv
grep -v -E 'NA|Inf|NaN' $OUTPATH/rrs_filter/TE_CDS_PC_teMerged.tsv >$OUTPATH/rrs_filter/TE_CDS_PC_teMerged_tidy.tsv

# Filter out rna_ratio_median < 45 and RRS_median < 5 (cutoff based on RRS paper)
#tail -n +2 RRS_out.tsv | awk '$4 != "Inf"' | awk '$4 > 5' >RRS_out_filtered.tsv
# NOTE: roughly 22% have very short UTR3/pUTR3, are they FP?
tail -n +2 $OUTPATH/rrs_filter/RRS_out.tsv | awk '$2 > 9 && $3 < 30 && $4 > 5' >$OUTPATH/rrs_filter/RRS_out_filtered.tsv

LL_ALL=C fgrep -f $OUTPATH/rrs_filter/ORFId_canonical.txt $OUTPATH/rrs_filter/RRS_out_filtered.tsv | cut -f1,4 | sort >$OUTPATH/rrs_filter/RRS_canonical.tsv
LL_ALL=C fgrep -v -f $OUTPATH/rrs_filter/ORFId_canonical.txt $OUTPATH/rrs_filter/RRS_out_filtered.tsv | cut -f1,4 | sort >$OUTPATH/rrs_filter/RRS_nonCanonical.tsv
cat $OUTPATH/rrs_filter/RRS_canonical.tsv $OUTPATH/rrs_filter/RRS_nonCanonical.tsv | sort -k2 -gr >$OUTPATH/rrs_filter/RRS_all.tsv
cut -f1 $OUTPATH/rrs_filter/RRS_all.tsv | sort >$OUTPATH/rrs_filter/ORFId_rrs_filter.txt

echo "Finished rrs filter."
echo

# Step 9:
# Nested filter, if there are ORFs with the same stop but different start, decide which to keep

echo "Start nested filter..."

mkdir -p $OUTPATH/nested_filter

LC_ALL=C fgrep -f $OUTPATH/rrs_filter/ORFId_rrs_filter.txt $OUTPATH/orfscore_filter/ORFScore_all_filteredByCov.tsv | sort -k1 >$OUTPATH/nested_filter/ORFScore_info.tsv
LC_ALL=C fgrep -f $OUTPATH/rrs_filter/ORFId_rrs_filter.txt $OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed | sort -k4 >$OUTPATH/nested_filter/BED_info.tsv
paste $OUTPATH/nested_filter/BED_info.tsv $OUTPATH/nested_filter/ORFScore_info.tsv | cut -f1-13,15-26 >$OUTPATH/nested_filter/BED12Plus.bed

Rscript $CODEBASE/util/NestedRegionFilter.R $OUTPATH/nested_filter $THREAD

cut -f6 $OUTPATH/nested_filter/nested_filtered.tsv | sort | uniq >$OUTPATH/nested_filter/ORFId_nested_filter.txt

echo "Finished nested filter."
echo

# Step 10:
# FDR filter, calulate FDR for ORFScore

echo "Start FDR filter..."

mkdir -p $OUTPATH/fdr_filter

LC_ALL=C fgrep -f $OUTPATH/nested_filter/ORFId_nested_filter.txt $OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed | cut -f4,13 | sort >$OUTPATH/fdr_filter/ORFId_label.tsv
LC_ALL=C fgrep -f $OUTPATH/nested_filter/ORFId_nested_filter.txt $OUTPATH/orfscore_filter/ORFScore_all_noNA.tsv | sort >$OUTPATH/fdr_filter/ORFId_orfscore.tsv 
LC_ALL=C fgrep -f $OUTPATH/nested_filter/ORFId_nested_filter.txt $OUTPATH/region_filter/orfs_readRatio.tsv | cut -f6,15 | sort >$OUTPATH/fdr_filter/ORFId_regionFilter.tsv
LC_ALL=C fgrep -f $OUTPATH/nested_filter/ORFId_nested_filter.txt $OUTPATH/rrs_filter/RRS_all.tsv | sort >$OUTPATH/fdr_filter/RRS_fdr.tsv
paste $OUTPATH/fdr_filter/ORFId_orfscore.tsv $OUTPATH/fdr_filter/ORFId_label.tsv $OUTPATH/fdr_filter/RRS_fdr.tsv $OUTPATH/fdr_filter/ORFId_regionFilter.tsv | cut -f1-13,15,17,19 | sort -k2 -g >$OUTPATH/fdr_filter/ORFScore_merged.tsv

# Remove retained_intron
grep -v $OUTPATH/fdr_filter/retained_intron $OUTPATH/fdr_filter/ORFScore_merged.tsv >$OUTPATH/fdr_filter/ORFScore_merged_RI_filtered.tsv

# Filter for processed_transcript, FP cases are those overlapping with pc transcripts
grep processed_transcript $OUTPATH/fdr_filter/ORFScore_merged_RI_filtered.tsv | awk '$16 != "NA"' | cut -f1 >$OUTPATH/fdr_filter/ORFId_processed_transcript_filter_1.txt
grep processed_transcript $OUTPATH/fdr_filter/ORFScore_merged_RI_filtered.tsv | awk '$14 != "noncoding"' | cut -f1 >$OUTPATH/fdr_filter/ORFId_processed_transcript_filter_2.txt
cat $OUTPATH/fdr_filter/ORFId_processed_transcript_filter_1.txt $OUTPATH/fdr_filter/ORFId_processed_transcript_filter_2.txt | sort | uniq >$OUTPATH/fdr_filter/ORFId_processed_transcript_filter_1_2.txt
LC_ALL=C fgrep -v -f $OUTPATH/fdr_filter/ORFId_processed_transcript_filter_1_2.txt $OUTPATH/fdr_filter/ORFScore_merged_RI_filtered.tsv >$OUTPATH/fdr_filter/ORFScore_merged_PT_filtered_tmp.tsv
grep processed_transcript $OUTPATH/fdr_filter/ORFScore_merged_PT_filtered_tmp.tsv | cut -f1 >$OUTPATH/fdr_filter/ORFId_processed_transcript_tmp.txt
grep processed_transcript $OUTPATH/fdr_filter/ORFScore_merged_PT_filtered_tmp.tsv | cut -f2 -d: | sort | uniq >$OUTPATH/fdr_filter/ORFId_processed_transcript_geneSymbol_tmp.txt
LC_ALL=C fgrep -f $OUTPATH/fdr_filter/ORFId_processed_transcript_geneSymbol_tmp.txt $GTF | awk '$3 == "gene"' | cut -f9 | sed 's/"//g;s/;//g' | cut -f6,4 -d" " | awk '{print $2"\t"$1}' >$OUTPATH/fdr_filter/ORFId_processed_transcript_GeneSymbol_GeneType_tmp.tsv 
grep protein_coding $OUTPATH/fdr_filter/ORFId_processed_transcript_GeneSymbol_GeneType_tmp.tsv | cut -f1 | sort >$OUTPATH/fdr_filter/pc_symbol.txt
LC_ALL=C fgrep -f $OUTPATH/fdr_filter/pc_symbol.txt $OUTPATH/fdr_filter/ORFId_processed_transcript_tmp.txt >$OUTPATH/fdr_filter/ORFId_processed_transcript_pc_tmp.txt
LC_ALL=C fgrep -v -f $OUTPATH/fdr_filter/ORFId_processed_transcript_pc_tmp.txt $OUTPATH/fdr_filter/ORFScore_merged_PT_filtered_tmp.tsv >$OUTPATH/fdr_filter/ORFScore_merged_PT_filtered.tsv

# p-val adjust
Rscript $CODEBASE/util/ORFScore_padj.R $OUTPATH/fdr_filter ORFScore_merged_PT_filtered.tsv all_withQval

# q-val cutoff 0.01. If comparing numbers with E notation, use ($4+0), otherwise awk will think it is a string 
awk '($4+0) < 0.01' $OUTPATH/fdr_filter/ORFScore_all_withQval.tsv >$OUTPATH/fdr_filter/ORFScore_fdr_filtered.tsv

# Length filter >= 11aa including stop codon
cat $OUTPATH/fdr_filter/ORFScore_fdr_filtered.tsv | cut -f1 | cut -f4,5 -d: | tr ":" "\t" | awk '{print ($2-$1+1)/3}' >$OUTPATH/fdr_filter/ORFScore_fdr_filtered_len.txt
paste $OUTPATH/fdr_filter/ORFScore_fdr_filtered.tsv $OUTPATH/fdr_filter/ORFScore_fdr_filtered_len.txt | awk '$17 >= 11' | cut -f1-16 >$OUTPATH/fdr_filter/ORFScore_length_filtered.tsv

# codon_1_2 mean should be greater than 0
awk '$10 > 0' $OUTPATH/fdr_filter/ORFScore_length_filtered.tsv >$OUTPATH/fdr_filter/ORFScore_codon_filtered.tsv

cut -f1 $OUTPATH/fdr_filter/ORFScore_codon_filtered.tsv >$OUTPATH/fdr_filter/ORFId_fdr_filter.txt

echo "Finished FDR filter."
echo

# Step 11:
# Create output file

echo "Start creating output files..."

mkdir -p $OUTPATH/final_output

LC_ALL=C fgrep -f $OUTPATH/fdr_filter/ORFId_fdr_filter.txt $OUTPATH/label_filter/regions.tsv | sort -k2 >$OUTPATH/final_output/region_ORFId.tsv
LC_ALL=C fgrep -f $OUTPATH/fdr_filter/ORFId_fdr_filter.txt $OUTPATH/label_filter/regions.tsv | cut -f1 | sort | uniq >$OUTPATH/final_output/region_nonDup.txt

LC_ALL=C fgrep -f $OUTPATH/fdr_filter/ORFId_fdr_filter.txt $OUTPATH/label_filter/ORFScore_ALL_filtered_expressed_label_BED12Plus.bed | sort -k4 >$OUTPATH/final_output/smORFs_ALL_BED12Plus.bed
cut -f4 $OUTPATH/final_output/smORFs_ALL_BED12Plus.bed | cut -f2 -d: | sort | uniq >$OUTPATH/final_output/smORFs_gene_symbol.txt

paste $OUTPATH/final_output/smORFs_ALL_BED12Plus.bed $OUTPATH/final_output/region_ORFId.tsv | cut -f1-14 >$OUTPATH/final_output/smORFs_ALL_BED12Plus_withRegionId.bed

awk -v OFS="\t" '{print $1,$2,$3,$14,$5,$6,$7,$8,$9,$10,$11,$12,$13}' $OUTPATH/final_output/smORFs_ALL_BED12Plus_withRegionId.bed | sort | uniq >$OUTPATH/final_output/smORFs_ALL_BED12Plus_regionIdAsCol4.bed

LC_ALL=C fgrep -f $OUTPATH/fdr_filter/ORFId_fdr_filter.txt $OUTPATH/fdr_filter/ORFScore_merged.tsv | cut -f14 | sort | uniq | grep -v '^$' >$OUTPATH/final_output/label.txt
parallel -j $THREAD "grep -P \"{}$\" $OUTPATH/final_output/smORFs_ALL_BED12Plus_regionIdAsCol4.bed > $OUTPATH/final_output/smORFs_{}_BED12Plus.bed" :::: $OUTPATH/final_output/label.txt

parallel -j $THREAD "cut -f1-12 $OUTPATH/final_output/smORFs_{}_BED12Plus.bed >$OUTPATH/final_output/smORFs_{}_BED12.bed" :::: $OUTPATH/final_output/label.txt

parallel -j $THREAD "bedtools getfasta -fi $REFGENOME -bed $OUTPATH/final_output/smORFs_{}_BED12.bed -split -name -s -fo $OUTPATH/final_output/smORFs_{}.fa" :::: $OUTPATH/final_output/label.txt

parallel -j $THREAD "awk '/^>/{print \$1\" \"\$1; next}{print}' <{} >{.}_reHeader.fa" ::: `ls $OUTPATH/final_output/*.fa`

parallel -j $THREAD "transeq {} {.}.pep.t -frame=1" ::: `ls $OUTPATH/final_output/*_reHeader.fa`

parallel -j $THREAD "awk '/^>/{print \$2; next}{print}' <{} >{.}" ::: `ls $OUTPATH/final_output/*pep.t`
rm $OUTPATH/final_output/*.t
rm $OUTPATH/final_output/smORFs_*_reHeader.fa

cat $OUTPATH/final_output/*.pep | sed 's/(.)//g' >$OUTPATH/final_output/all.pep
cat $OUTPATH/final_output/all.pep | awk '$0 ~ ">" {print c; c=0;printf substr($0,2) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sed '/^$/d' | awk '{print $1"\t"$2-1}' | sort | uniq >$OUTPATH/final_output/pep_len.tsv 

echo "Finished creating output files."
echo

# Step 12:
# Create information table

echo "Start creating information table..."

mkdir -p $OUTPATH/info_table

cut -f4 $OUTPATH/final_output/smORFs_ALL_BED12Plus.bed | sort >$OUTPATH/info_table/ORFId_final.txt
cut -f1 -d: $OUTPATH/info_table/ORFId_final.txt >$OUTPATH/info_table/host_tx.txt
cut -f2 -d: $OUTPATH/info_table/ORFId_final.txt >$OUTPATH/info_table/host_gene.txt

paste $OUTPATH/info_table/host_tx.txt $OUTPATH/info_table/ORFId_final.txt $OUTPATH/info_table/host_gene.txt | sort >$OUTPATH/info_table/Ids.tsv

# extract ORFscore and RRS
LC_ALL=C fgrep -f $OUTPATH/info_table/ORFId_final.txt $OUTPATH/fdr_filter/ORFScore_codon_filtered.tsv | cut -f1-2,15 | sort >$OUTPATH/info_table/ORFId_ORFScore_RRS.tsv

# extract RPKM, TE_smORF and TE_mainORF
LC_ALL=C fgrep -f $OUTPATH/info_table/ORFId_final.txt $OUTPATH/rrs_filter/TE_smORF_full.tsv | cut -f2,5,6 | sort >$OUTPATH/info_table/ORFId_RPrpkm_RNArpkm.tsv

LC_ALL=C fgrep -f $OUTPATH/info_table/ORFId_final.txt $OUTPATH/rrs_filter/TE_smORF_full.tsv | cut -f2,7 | sort >$OUTPATH/info_table/ORFId_TEsmORF.tsv

# Extract both uORF and ouORF
grep five_prime $OUTPATH/final_output/smORFs_ALL_BED12Plus.bed | cut -f4 >$OUTPATH/info_table/five_prime_ORFId.txt
cat $OUTPATH/info_table/five_prime_ORFId.txt | cut -f1 -d: >$OUTPATH/info_table/five_prime_host_txId.txt
paste $OUTPATH/info_table/five_prime_ORFId.txt $OUTPATH/info_table/five_prime_host_txId.txt >$OUTPATH/info_table/five_prime_ORFId_txId.tsv 
cat $OUTPATH/info_table/five_prime_host_txId.txt | sort | uniq >$OUTPATH/info_table/five_prime_host_txId_uniq.txt


LC_ALL=C fgrep -f $OUTPATH/info_table/five_prime_host_txId_uniq.txt $OUTPATH/rrs_filter/TE_CDS_PC_teMerged.tsv | cut -f1,2 | sort >$OUTPATH/info_table/TxId_TEmainORF.tsv


join -1 2 -2 1 <(sort -k2 $OUTPATH/info_table/five_prime_ORFId_txId.tsv) $OUTPATH/info_table/TxId_TEmainORF.tsv | awk -v OFS='\t' '{print $2,$3}' | sort >$OUTPATH/info_table/five_prime_ORFId_TEmainORF.tsv

paste $OUTPATH/info_table/ORFId_ORFScore_RRS.tsv $OUTPATH/info_table/ORFId_RPrpkm_RNArpkm.tsv $OUTPATH/info_table/ORFId_TEsmORF.tsv | cut -f1-3,5,6,8 >$OUTPATH/info_table/ORFId_ORFScore_RRS_RPrpkm_RNArpkm_TEsmORF.tsv

join -t $'\t' -a1 -a2 -e 'NA' -o '0,1.2,1.3,1.4,1.5,1.6,2.2' $OUTPATH/info_table/ORFId_ORFScore_RRS_RPrpkm_RNArpkm_TEsmORF.tsv $OUTPATH/info_table/five_prime_ORFId_TEmainORF.tsv | sort >$OUTPATH/info_table/ORFId_ORFScore_RRS_RPrpkm_RNArpkm_TEsmORF_TEmainORF.tsv

# extract gene description
sort $OUTPATH/info_table/host_tx.txt | uniq >$OUTPATH/info_table/host_tx_uniq.txt 
LC_ALL=C fgrep -f $OUTPATH/info_table/host_tx_uniq.txt $GENEDESCRIPTION | sort >$OUTPATH/info_table/host_tx_idsWithVersion_geneName_description.tsv

join -t $'\t' -a1 -a2 $OUTPATH/info_table/Ids.tsv $OUTPATH/info_table/host_tx_idsWithVersion_geneName_description.tsv | awk -F '\t' -v OFS="\t" '{if ($3!=$4) print $2,$1,$3"|"$4,$5; else print $2,$1,$3,$5}' | sort >$OUTPATH/info_table/ORFId_txId_geneName_description.tsv

paste <(sort -k14 $OUTPATH/final_output/smORFs_ALL_BED12Plus_withRegionId.bed | cut -f14) <(sort -k14 $OUTPATH/final_output/smORFs_ALL_BED12Plus_withRegionId.bed | cut -f1-13) >$OUTPATH/info_table/smORF_withRegion_rearrange.tsv
join -t $'\t' -a1 -a2 $OUTPATH/info_table/smORF_withRegion_rearrange.tsv <(sort $OUTPATH/final_output/pep_len.tsv) >$OUTPATH/info_table/smORF_withRegion_pepLen.tsv
paste <(cut -f2-15 $OUTPATH/info_table/smORF_withRegion_pepLen.tsv) <(cut -f1 $OUTPATH/info_table/smORF_withRegion_pepLen.tsv) | sort -k4 >$OUTPATH/info_table/smORF_BED12_class_peplen_regionId.bed

paste <(cut -f4 $OUTPATH/info_table/smORF_BED12_class_peplen_regionId.bed) $OUTPATH/info_table/smORF_BED12_class_peplen_regionId.bed >$OUTPATH/info_table/smORF_BED12_class_peplen_regionId.tsv
join -t $'\t' -a1 -a2 $OUTPATH/info_table/smORF_BED12_class_peplen_regionId.tsv $OUTPATH/info_table/ORFId_txId_geneName_description.tsv >$OUTPATH/info_table/smORF_BED12_class_peplen_regionId_txId_geneName_description.tsv
join -t $'\t' -a1 -a2 $OUTPATH/info_table/smORF_BED12_class_peplen_regionId_txId_geneName_description.tsv $OUTPATH/info_table/ORFId_ORFScore_RRS_RPrpkm_RNArpkm_TEsmORF_TEmainORF.tsv | cut -f2-25 >$OUTPATH/info_table/info_table_BED12Plus_withoutAASeq.bed

perl -pe 's/>(.*)/>\1\t/g; s/\n//g; s/>/\n>/g' $OUTPATH/final_output/all.pep | grep -v '^$' | cut -c 2- | sort | uniq >$OUTPATH/info_table/RegionID_AASeq.tsv
paste <(cut -f15 $OUTPATH/info_table/info_table_BED12Plus_withoutAASeq.bed) $OUTPATH/info_table/info_table_BED12Plus_withoutAASeq.bed | sort >$OUTPATH/info_table/RegionId_info_table_BED12Plus_withoutAASeq.tsv
join -t $'\t' -a1 -a2 $OUTPATH/info_table/RegionId_info_table_BED12Plus_withoutAASeq.tsv $OUTPATH/info_table/RegionID_AASeq.tsv | cut -f2-26 | sort -n >$OUTPATH/info_table/info_table_BED12Plus.bed


echo "Finished creating information table."
echo

echo "Finished ORF calling."
echo
