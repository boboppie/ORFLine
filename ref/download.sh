# Download reference genome and annotations

##### GENCODE #####

VERSION=8
mkdir -p GENCODE/mouse/M$VERSION
GPATH=GENCODE/mouse/M$VERSION

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
