# Download reference genome and annotations

##### GENCODE #####

VERSION=8
mkdir -p ~/ref/GENCODE/mouse/M$VERSION
GPATH=~/ref/GENCODE/mouse/M$VERSION

### Sequence fasta
# Genome PRI (CHR + scaffolds)

mkdir -p $GPATH/fasta/genome/PRI
wget -O $GPATH/fasta/genome/PRI/GRCm38.primary_assembly.genome.fa.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/GRCm38.primary_assembly.genome.fa.gz

# Transcripts CHR

mkdir -p $GPATH/fasta/transcriptome/CHR/all-transcripts
mkdir -p $GPATH/fasta/transcriptome/CHR/protein-coding-transcripts
mkdir -p $GPATH/fasta/transcriptome/CHR/lncRNA-transcripts

wget -O $GPATH/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.transcripts.fa.gz
wget -O $GPATH/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.pc_transcripts.fa.gz
wget -O $GPATH/fasta/transcriptome/CHR/lncRNA-transcripts/lncRNA_transcripts.fa.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.lncRNA_transcripts.fa.gz

### Annitation GTF/GFF3
mkdir -p $GPATH/annotation/CHR/comprehensive
mkdir -p $GPATH/annotation/CHR/basic
mkdir -p $GPATH/annotation/CHR/lncRNA
mkdir -p $GPATH/annotation/CHR/polyAs
mkdir -p $GPATH/annotation/CHR/pseudogenes
mkdir -p $GPATH/annotation/CHR/tRNA

# Comprehensive gene annotation CHR

wget -O $GPATH/annotation/CHR/comprehensive/annotation.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.annotation.gtf.gz
wget -O $GPATH/annotation/CHR/comprehensive/annotation.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.annotation.gff3.gz

# Basic gene annotation CHR

wget -O $GPATH/annotation/CHR/basic/basic.annotation.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.basic.annotation.gtf.gz
wget -O $GPATH/annotation/CHR/basic/basic.annotation.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.basic.annotation.gff3.gz

# Long non-coding RNA gene annotation CHR

wget -O $GPATH/annotation/CHR/lncRNA/lncRNA.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.long_noncoding_RNAs.gtf.gz
wget -O $GPATH/annotation/CHR/lncRNA/lncRNA.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.long_noncoding_RNAs.gff3.gz

# PolyA feature annotation CHR

wget -O $GPATH/annotation/CHR/polyAs/polyAs.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.polyAs.gtf.gz
wget -O $GPATH/annotation/CHR/polyAs/polyAs.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.polyAs.gff3.gz

# Consensus pseudogenes predicted by the Yale and UCSC pipelines CHR

wget -O $GPATH/annotation/CHR/pseudogenes/pseudos.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.2wayconspseudos.gtf.gz
wget -O $GPATH/annotation/CHR/pseudogenes/pseudos.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.2wayconspseudos.gff3.gz

# Predicted tRNA genes CHR 

wget -O $GPATH/annotation/CHR/tRNA/tRNAs.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.tRNAs.gtf.gz
wget -O $GPATH/annotation/CHR/tRNA/tRNAs.gff3.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.tRNAs.gff3.gz

### Metadata ALL

# Entrez gene id

mkdir -p $GPATH/metadata/ALL
wget -O $GPATH/metadata/ALL/EntrezGene.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.metadata.EntrezGene.gz

# Gene symbol

wget -O $GPATH/metadata/ALL/MGI.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.metadata.MGI.gz

# RefSeq

wget -O $GPATH/metadata/ALL/RefSeq.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M$VERSION/gencode.vM$VERSION.metadata.RefSeq.gz

###########################
gunzip $GPATH/fasta/genome/PRI/GRCm38.primary_assembly.genome.fa.gz
mkdir -p $GPATH/fasta/genome/CHR
~/code/github/orf-discovery/script/fastagrep.sh 'chr' $GPATH/fasta/genome/PRI/GRCm38.primary_assembly.genome.fa > $GPATH/fasta/genome/CHR/GRCm38.CHR.genome.fa
bowtie-build $GPATH/fasta/genome/CHR/GRCm38.CHR.genome.fa $GPATH/fasta/genome/CHR/GRCm38.CHR.genome
bowtie2-build $GPATH/fasta/genome/CHR/GRCm38.CHR.genome.fa $GPATH/fasta/genome/CHR/GRCm38.CHR.genome
samtools faidx $GPATH/fasta/genome/CHR/GRCm38.CHR.genome.fa

# RTCm38/mm10 tRNAs 
# GENCODE vM8 has 26248 entries predicted by tRNAscan-SE
# GtRNAdb 432 (ref http://lowelab.ucsc.edu/GtRNAdb/Mmusc10/)
# UCSC 435(ref http://onetipperday.sterding.com/2012/08/how-to-get-trnarrnamitochondrial-gene.html)
# group: All Tables
# database: mm10
# table: tRNAs
# Download (manually) from UCSC Table Browser UCSC.tRNAs.fa to ~/ref/UCSC/mouse/mm10/fasta/tRNA/
# Remove (manually) the only chrM entry in UCSC >mm10_tRNAs_chrM.tRNA1-LeuTAA range=chrM:2676-2750 which is mt-Tl1 in GENCODE

# fasta for tRNA, rRNA, Mt_rRNA, Mt_tRNA, snRNA, snoRNA, misc_RNA, miRNA for ribo-seq reads filtering
mkdir -p $GPATH/fasta/ribo-seq_reads-filter
~/code/github/orf-discovery/script/fastagrep.sh 'rRNA|Mt_rRNA|Mt_tRNA|snRNA|snoRNA|misc_RNA|miRNA' <(gzip -dc $GPATH/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz) > $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter_notRNAs.fa # 6036 entries
cat $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter_notRNAs.fa ~/ref/UCSC/mouse/mm10/fasta/tRNA/UCSC.tRNAs.fa >$GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter.fa # 6470 in total
bowtie-build $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter.fa $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter
bowtie2-build $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter.fa $GPATH/fasta/ribo-seq_reads-filter/ribo-seq_reads-filter

mkdir -p $GPATH/fasta/ribo-seq_reads-unfiltered
~/code/github/orf-discovery/script/fastagrep.sh -v 'rRNA|Mt_rRNA|Mt_tRNA|snRNA|snoRNA|misc_RNA|miRNA' <(gzip -dc $GPATH/fasta/transcriptome/CHR/all-transcripts/transcripts.fa.gz) > $GPATH/fasta/ribo-seq_reads-unfiltered/ribo-seq_reads-unfiltered.fa
bowtie-build $GPATH/fasta/ribo-seq_reads-unfiltered/ribo-seq_reads-unfiltered.fa $GPATH/fasta/ribo-seq_reads-unfiltered/ribo-seq_reads-unfiltered
bowtie2-build $GPATH/fasta/ribo-seq_reads-unfiltered/ribo-seq_reads-unfiltered.fa $GPATH/fasta/ribo-seq_reads-unfiltered/ribo-seq_reads-unfiltered

# Filter out annotation
mkdir -p $GPATH/annotation/CHR/filtered
egrep -iv 'snRNA|snoRNA|Mt_rRNA|Mt_tRNA|misc_RNA|miRNA|rRNA' <(gzip -dc $GPATH/annotation/CHR/comprehensive/annotation.gtf.gz) | gzip >$GPATH/annotation/CHR/filtered/filtered.gtf.gz

# Protein coding gtf 
# gene feature will be filtered out, ref - https://www.biostars.org/p/101595/
# 50620 transcript entries, however in fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.gz, there are 56504 entries, e.g. ENSMUST00000006138.12 is nonsense_mediated_decay
#mkdir -p $GPATH/annotation/CHR/protein_coding
#zgrep 'transcript_type "protein_coding"' GPATH/annotation/CHR/comprehensive/annotation.gtf.gz | gzip >GPATH/annotation/CHR/protein_coding/protein_coding.gtf.gz

#gunzip ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.gz
#grep '>' ~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa | cut -f1 -d"|" | cut -f2 -d'>' | sort >~/ref/GENCODE/mouse/M8/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.tlst.ids

# In order to match protein_coding fasta, filtering based on GENCODE description: Transcript biotypes: protein_coding, nonsense_mediated_decay, non_stop_decay, IG_*_gene, TR_*_gene, polymorphic_pseudogene
egrep -i 'transcript_type "nonsense_mediated_decay"|transcript_type "protein_coding"|transcript_type "IG_.*_gene"|transcript_type "TR_.*_gene"|transcript_type "polymorphic_pseudogene"|transcript_type "non_stop_decay"' $GPATH/annotation/CHR/comprehensive/annotation.gtf >$GPATH/annotation/CHR/protein_coding/protein_coding.gtf 

awk '{if($3=="transcript"){print $0}}' $GPATH/annotation/CHR/protein_coding/protein_coding.gtf | cut -f9 | cut -f2 -d";" | cut -f3 -d" " | sed 's/"//g' | sort | uniq >$GPATH/annotation/CHR/protein_coding/protein_coding.gtf.tlst.id

diff $GPATH/annotation/CHR/protein_coding/protein_coding.gtf.tlst.id $GPATH/fasta/transcriptome/CHR/protein-coding-transcripts/pc_transcripts.fa.tlst.ids | grep '<' | cut -f2 -d" " >$GPATH/annotation/CHR/protein_coding/pc_gtf_fa.diff.ids

egrep -v `cat $GPATH/annotation/CHR/protein_coding/pc_gtf_fa.diff.ids | tr '\n' '|' | sed s'/.$//'` $GPATH/annotation/CHR/protein_coding/protein_coding.gtf >$GPATH/annotation/CHR/protein_coding/protein_coding.filtered.gtf

# Selenocysteine will break TopHat
grep Selenocysteine $GPATH/annotation/CHR/protein_coding/protein_coding.filtered.gtf >$GPATH/annotation/CHR/protein_coding/protein_coding.filtered.Selenocysteine.gtf
grep -v Selenocysteine $GPATH/annotation/CHR/protein_coding/protein_coding.filtered.gtf >$GPATH/annotation/CHR/protein_coding/protein_coding.filtered.Selenocysteine.gtf
