# Discovery of bioactive micropeptides in the immune system

This repository holds the pipeline for prediction of actively translated small open reading frames (smORFs) in the immune system.

## Obtaining

To download the source code, please use git to download the most recent development
tree.  Currently, the tree is hosted on github, and can be obtained via:

    git clone git://github.com/boboppie/orf-discovery.git
    
## Usage

### Dependencies

* [Samtools](http://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [BEDOPS](https://bedops.readthedocs.io/en/latest/)
* [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [STAR](https://github.com/alexdobin/STAR)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [plastid](https://plastid.readthedocs.io/en/latest/index.html)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [EMBOSS](http://emboss.sourceforge.net/)
* [GNU Parallel](https://www.gnu.org/software/parallel/) (recommended)
* [R](https://www.r-project.org/)

R/Bioconductor packages:
* [riboSeqR](http://bioconductor.org/packages/release/bioc/html/riboSeqR.html)
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

### Workflow

1. Download and generate files that are used in the pipeline

```bash
RUN_THE_SCRIPT
```

2. Generate putative ORFs

```bash
RUN_THE_SCRIPT
```

3. Ribosome profiling (Ribo-Seq) data processing

```bash
RUN_THE_SCRIPT
```

4. RNA-Seq data processing

```bash
RUN_THE_SCRIPT
```

5. ORF calling

```bash
RUN_THE_SCRIPT
```

## Pipeline explained

### Annotation

We recommend to use [GENCODE](https://www.gencodegenes.org/) annotation (GTF/GFF3 files) and genome and transcript sequences (Fasta files) for Human and Mouse. For other species, we have used [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) annotation for Zebrafish. tRNA sequences can be downloaded from UCSC Table Browser.

For example, we have downloaded the following public sequences and annotation files for Mouse:

File | Type | Region | Source
---- | ---- | ------ | ------
Genome sequence, primary assembly | Nucleotide sequence of the GRCm38 primary genome assembly (Fasta format) | PRI (reference chromosomes and scaffolds) | GENCODE
Transcript sequences | Nucleotide sequences of all transcripts (Fasta format) | CHR (reference chromosomes only) | GENCODE
Protein-coding transcript sequences | Nucleotide sequences of coding transcripts (Fasta format) | CHR | GENCODE
LncRNA transcript sequences | Nucleotide sequences of lncRNA transcripts (Fasta format) | CHR | GENCODE
Comprehensive gene annotation | The main annotation file (GTF and GFF3 format) | CHR | GENCODE
LncRNA gene annotation | comprehensive gene annotation of lncRNA genes (GTF and GFF3 format) | CHR | GENCODE
tRNA sequences | Nucleotide sequences of tRNA genes predicted by UCSC using tRNAscan-SE (Fasta format) | CHR | UCSC Table Browser

### Defining transcript sets

We combined GENCODE protein-coding transcripts and LncRNA transcripts to form reference transcriptome. Users can potentially assemble the transcriptome using RNA-seq data (e.g. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)), however, Studies have shown that computational approaches produce a large number of artefacts (false positives), which absorbed a substantial proportion of the reads from truly expressed transcripts and were assigned large expression estimates.

We do not consider the following [biotypes](https://www.gencodegenes.org/pages/biotypes.html):
* IG_* and TR_* (Immunoglobulin variable chain and T-cell receptor genes)
* miRNA
* misc_RNA
* Mt_rRNA and Mt_tRNA
* rRNA and ribozyme
* scaRNA, scRNA, snoRNA, snRNA and sRNA
* nonsense_mediated_decay
* non_stop_decay


For example, transcriptome can be generated by:

```bash
# Assume the source code is in ~/code/github/orf-discovery/
# Mouse GENCODE M13 reference files are in ~/ref/GENCODE/mouse/M13/
# Transcript sequence file downloaded from GENCODE is named gencode.vM13.transcripts.fa.gz, and placed in ~/ref/GENCODE/mouse/M13/fasta/transcriptome/
# The newly generated transcript sequence file name is transcripts_biotype_filtered.fa

cd ~/ref/GENCODE/mouse/M13/fasta/transcriptome/

~/code/github/orf-discovery/script/fastagrep.sh -v 'IG_*|TR_*|miRNA|misc_RNA|Mt_*|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|ribozyme|nonsense_mediated_decay|non_stop_decay' <(gzip -dc gencode.vM13.transcripts.fa.gz) >transcripts_biotype_filtered.fa
```

### Defining ORFs

Given transcriptome sequences, we exhaustively searched for putative ORFs beginning with a start codon (“ATG”, “TTG”, “CTG”, “GTG”) and ending with a stop codon ("TAG", "TAA", "TGA") without an intervening stop codon in between in each of the three reading frames.

Code to generate putative ORFs:

```bash
# Arguments:
# START_CODON - start codon, can be “ATG”, “TTG”, “CTG”, “GTG”  
# TRANSCRIPTOME_FASTA - transcriptome file location, assume is it in ~/ref/GENCODE/mouse/M13/fasta/transcriptome
# GTF - gene annotation file, assume it is in ~ref/GENCODE/mouse/M13/annotation/CHR/comprehensive/, and named gencode.vM13.annotation.gtf
# ORGANISM - scientific name, e.g. Mus musculus
# NCORE - number of computing cores if you want the process to be multi-threading, default 1
#
# for example, they can be assigned as: 
# START_CODON=ATG
# TRANSCRIPTOME_FASTA=~/ref/GENCODE/mouse/M13/fasta/transcriptome/transcripts_biotype_filtered.fa
# GTF=~ref/GENCODE/mouse/M13/annotation/CHR/comprehensive/gencode.vM13.annotation.gtf
# ORGANISM="Mus musculus"
# NCORE=8 
#
# Output:
# The output file is in BED12 format (https://genome.ucsc.edu/FAQ/FAQformat#format1) and named orf_$START_CODON.bed 

Rscript ~/code/github/orf-discovery/script/ORFPredict.R $START_CODON $TRANSCRIPTOME_FASTA $GTF $ORGANISM $NCORE 
```

In the output file, a unique ID is given for each ORF, for example:

    ENSMUST00000000010.8:Hoxb9:protein_coding:117:134:2:ATG

the fields separated by colons are **transcript ID**, **gene symbol**, **transcript biotype**, **transcript start**, **transcript stop**, **reading frame**, **start codon**. 

User can ran the R script for each start codon and keep them in the same directory.

We are interesed in smORFs (less than 100 codons). In order to filter for them, firstly, we calculate the length of ORFs, for example:

```bash
cut -f11 orfs_ATG.bed | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1' >orfs_ATG.width
```

Then, we only keep the smORFs:

```bash
paste -d'\t' orfs_ATG.bed orfs_ATG.width | awk '$13 <= 303' | cut -f1-12 | sort -k4 >orfs_ATG.smORFs_100.bed
```

### Ribo-Seq data processing

Ribo-Seq libraries are in general single-end and reads are 50bp long. We have the following steps to process the raw sequencing data (Fastq format):
    
1. Quality control. We use FastQC, for example: 

```bash
fastqc -t $THREADS -o $OUTPATH/fastqc $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz
```

2. Adapter and quality trimming. We use Trim Galore ([user guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)), for example:

```bash
# User can give a customized adapter by flag -a, or trim_galore will automatically detect the adapter if it's standard
# To filter for read length, user can use flags such as --length 25 --max_length 35 for min and max length

trim_galore -q 33 --fastqc --trim-n -e 0.1 --stringency 3 $RAWDATAPATH/${RAWDATAFILENAME}.fastq.gz -o $OUTPATH/trim_galore
```

3. Contaimination removing

   In order to remove rRNA/tRNA content or other contaminants in the sample, we used Bowtie (version 1) to align the trimmed  reads against specific contaminant sequences assembled from a collection of rRNA, Mt_rRNA, Mt_tRNA, snRNA, snoRNA, misc_RNA, miRNA (from GENCODE) and tRNA (from UCSC) sequences, we also include the following sequences from NCBI:  
    
       TPA_exp: Mus musculus ribosomal DNA, complete repeating unit
       http://www.ncbi.nlm.nih.gov/nuccore/511668571?report=fasta

       Mus musculus 45S pre-ribosomal RNA (Rn45s), ribosomal RNA
       http://www.ncbi.nlm.nih.gov/nuccore/577019615?report=fasta

       Mus musculus 28S ribosomal RNA (Rn28s1), ribosomal RNA
       http://www.ncbi.nlm.nih.gov/nuccore/120444900?report=fasta

       Mus musculus strain BALB/c 45S ribosomal RNA region genomic sequence
       http://www.ncbi.nlm.nih.gov/nuccore/307829144?report=fasta

       Mus musculus 4.5s RNA, pseudogene 1 (Rn4.5s-ps1) on chromosome 1
       http://www.ncbi.nlm.nih.gov/nuccore/693074770?report=fasta

       Microarray spike-in control plasmid pNIAysic-5, complete sequence
       http://www.ncbi.nlm.nih.gov/nuccore/70672673?report=fasta

       Human 28S ribosomal RNA gene
       http://www.ncbi.nlm.nih.gov/nuccore/337381?report=fasta

       Human 28S ribosomal RNA gene, complete cds
       http://www.ncbi.nlm.nih.gov/nuccore/337384?report=fasta

       Human ribosomal DNA complete repeating unit
       http://www.ncbi.nlm.nih.gov/nuccore/555853?report=fasta

       Homo sapiens RNA, 45S pre-ribosomal 5 (RNA45S5), ribosomal RNA
       http://www.ncbi.nlm.nih.gov/nuccore/374429547?report=fasta

       Chain 5, Structure Of The H. Sapiens 60s Rrna
       http://www.ncbi.nlm.nih.gov/nuccore/485601478?report=fasta

       NR_046235.3 Homo sapiens RNA, 45S pre-ribosomal (LOC100861532), ribosomal RNA
       NR_137294.1 Homo sapiens mitochondrially encoded 12S ribosomal RNA (RNR1), ribosomal RNA
       NR_137295.1 Homo sapiens mitochondrially encoded 16S ribosomal RNA (RNR2), ribosomal RNA
    
   Then Bowtie index will be built for the sequences:
   
   ```bash
   bowtie 
   ```
    
4. align to reference genome
5. P-site calling

### RNA-Seq data processing

QC
adapter/quality trimming
align to reference genome

### ORF calling

## Support

### email

Please report any issues or questions by creating a ticket, or by email to 
<fengyuan.hu@babraham.ac.uk>.
