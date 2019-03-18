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
* [R](https://www.r-project.org/)

R/Bioconductor packages:
* [riboSeqR](http://bioconductor.org/packages/release/bioc/html/riboSeqR.html)
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

### Tutorial to run the pipeline

#### Annotation

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

#### Defining transcript sets

We combined GENCODE protein-coding transcripts and LncRNA transcripts to form reference transcriptome. Users can potentially assemble the transcriptome using RNA-seq data (e.g. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)), however, Studies have shown that computational approaches produce a large number of artefacts (false positives), which absorbed a substantial proportion of the reads from truly expressed transcripts and were assigned large expression estimates.

#### Defining ORFs

Given transcriptome sequences, we exhaustively searched for ORFs beginning with a start codon (“ATG”, “TTG”, “CTG”, “GTG”) and ending with a stop codon ("TAG", "TAA", "TGA") without an intervening stop codon in between in each of the three reading frames.

## Support

### email

Please report any issues or questions by creating a ticket, or by email to 
<fengyuan.hu@babraham.ac.uk>.
