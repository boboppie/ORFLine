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

### Tutorial to run the pipeline

#### Annotation

We recommend to use [GENCODE](https://www.gencodegenes.org/) annotation (GTF/GFF3 files) and genome and transcript sequences (FASTA files) for Human and Mouse. For other species, we have used [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) annotation for Zebrafish.



## Support

### email

Please report any issues or questions by creating a ticket, or by email to 
<fengyuan.hu@babraham.ac.uk>.
