# Prediction of bioactive micropeptides in the immune system

This repository holds the pipeline for prediction of actively translated small open reading frames (smORFs) in the immune system.

## Obtaining

To download the source code, please use git to download the most recent development
tree.  Currently, the tree is hosted on github, and can be obtained via:

    git clone git://github.com/boboppie/orf-discovery.git
    
## Usage

### Dependencies

* [Samtools and HTSlib](http://www.htslib.org/)
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
* [Bioconductor](http://bioconductor.org/)

R/Bioconductor packages:
* [riboSeqR](http://bioconductor.org/packages/release/bioc/html/riboSeqR.html)
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

### Dataset

We will use [*Diaz-Muñoz et al, 2015*](https://www.nature.com/articles/ni.3115) LPS activated B cell dataset as an example to demonstrate typical workflow. 

Download raw sequencing data from EBI:

    RNA-Seq  -  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR160/001/SRR1605271/SRR1605271.fastq.gz
    Ribo-Seq -  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR160/004/SRR1605304/SRR1605304.fastq.gz

### Workflow

1. Check if all the dependencies are isntalled

```bash
bash ./module-check.sh
```

2. Download and generate files that are used in the pipeline

```bash
bash ./ref-download.sh -o mouse -r M22 -t 4
```

3. Generate putative ORFs

```bash
bash ./orf-prediction.sh -o \"Mus musculus\" -t 8
```

4. Ribosome profiling (Ribo-Seq) data processing

```bash
bash ./riboseq-process.sh -f ./out/data/ribo-seq/ribo.fastq.gz -a AAAAAAAAAAA -t 4
```

5. RNA-Seq data processing

```bash
bash ./rnaseq-process.sh -f ./out/data/rna-seq/rna.fastq.gz -t 4
```

6. ORF calling

```bash
bash ./orf-calling.sh -o mouse -x 10090 -m 32 -n 28 -t 8
```

### Output

The final output file in info_table directory is in BED12 format with extension.

column | Description
------ | -----------
1 - 12 | See [BED format definition](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). Col 4 is the ORFId
13 | smORF class, e.g. canonical, five_prime.
14 | Peptide length (AA)
15 | RegionId, it is possible muliple ORFIds (transcript-based) map to a unique regionId (genomic-based)
16 | Ensembl transcript Id
17 | Gene symbol
18 | Gene description
19 | ORF score
20 | Ribosome release score
21 | Ribo FPFM
22 | RNA FPKM
23 | Translation efficiency (TE)
24 | Main ORF TE (NA if host transcript is noncoding)


### Run the pipeline

We recommand to run a test on a virtual machine, e.g. [VirtualBox](https://www.virtualbox.org/wiki/Downloads) with a minimal ISO (e.g. [CentOS 7 minimal](http://isoredirect.centos.org/centos/7/isos/x86_64/CentOS-7-x86_64-Minimal-1810.iso)). Users can install all dependencies via [miniconda](https://docs.conda.io/en/latest/miniconda.html), for example:

```bash
# Tools to install on CentOS before miniconda
yum -y install gcc tar bzip2 git which

curl -fsSL https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -o miniconda2.sh

# assume miniconda is installed in the home directory
bash miniconda2.sh -b -p ~/miniconda2

export PATH=~/miniconda2/bin:$PATH
export PYTHONPATH=~/miniconda2/lib/python2.7/site-packages

conda install -y -c conda-forge wget 
conda install -y -c conda-forge parallel
conda install -y -c bioconda samtools
conda install -y -c bioconda htslib 
conda install -y -c bioconda bedtools 
conda install -y -c bioconda bedops 
conda install -y -c bioconda bowtie 
conda install -y -c bioconda fastqc
conda install -y -c bioconda cutadapt 
conda install -y -c bioconda trim-galore 
conda install -y -c bioconda star 
conda install -y -c bioconda stringtie 
conda install -y -c bioconda sra-tools
conda install -y -c bioconda emboss 
conda install -y -c bioconda plastid
conda install -y -c bioconda bioconductor-rhtslib
Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org"); BiocManager::install(c("riboSeqR", "GenomicFeatures", "rtracklayer"))'
```

We have a main.sh script to run all the steps mentioned above, you can simply pull the source code and run it as:

```bash
git clone https://github.com/boboppie/orf-discovery.git
cd orf-discovery
chmod +x *.sh
  
bash ./main.sh
```

### Singularity Container

We have created a [Singularity](https://singularity-hub.org) image, the repository is [here](https://github.com/boboppie/orf-discovery-singularity), users can pull and run a test as:

```bash
singularity pull shub://boboppie/orf-discovery-singularity # the default name of the image is orf-discovery-singularity_latest.sif

singularity run orf-discovery-singularity_latest.sif
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

Example options to download mouse tRNA sequences from [UCSC table brower](https://genome.ucsc.edu/cgi-bin/hgTables):
* clade:  Mammal     
* genome: Mouse
* assembly:   Dec 2011. (GRCm38/mm10)
* group:  All Tracks
* track:  tRNA Genes
* table:  tRNAs
* output format:  sequence

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

~/code/github/orf-discovery/script/fastagrep.sh -v 'IG_*|TR_*|miRNA|misc_RNA|Mt_*|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|ribozyme|nonsense_mediated_decay|non_stop_decay' <(gzip -dc <gencode_transcripts.fa.gz>) > <transcriptome.fa>
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

Rscript ~/code/github/orf-discovery/script/ORFPredict.R <start_codon> <transcriptome.fa> <annotation.gtf> <organism_scientific_name> <threads> 
```

In the output file, a unique ID is given for each ORF, for example:

    ENSMUST00000000010.8:Hoxb9:protein_coding:117:134:2:ATG

the fields separated by colons are **transcript ID**, **gene symbol**, **transcript biotype**, **transcript start**, **transcript stop**, **reading frame**, **start codon**. 

User can ran the R script for each start codon and keep them in the same directory.

We are interesed in smORFs (less than 100 codons). In order to filter for them, firstly, we calculate the length of ORFs, for example:

```bash
cut -f11 <orfs_ATG.bed> | sed -e 's/,/\t/g' | awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1' > <orfs_ATG.width>
```

Then, we only keep the smORFs:

```bash
paste -d'\t' <orfs_ATG.bed> <orfs_ATG.width> | awk '$13 <= 303' | cut -f1-12 | sort -k4 > <orfs_ATG.smORFs_100.bed>
```

### Ribo-Seq data processing

Ribo-Seq libraries are in general single-end and reads are 50bp long. We have the following steps to process the raw sequencing data (Fastq format):
    
#### 1. Quality control. We use FastQC, for example: 

   ```bash
   fastqc -t <threads> -o <fastqc_dir> <sample.fastq.gz>
   ```

#### 2. Adapter and quality trimming. We use Trim Galore ([user guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)), for example:

   ```bash
   # User can give a customized adapter by flag -a, or trim_galore will automatically detect the adapter if it's standard
   # To filter for read length, user can use flags such as --length 25 --max_length 35 for min and max length

   trim_galore -q 33 --fastqc --trim-n -e 0.1 --stringency 3 <sample.fastq.gz> -o <trim_galore_dir>
   ```

#### 3. Contaimination removing

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
       https://www.ncbi.nlm.nih.gov/nuccore/NR_046235.3/?report=fasta
       
       NR_137294.1 Homo sapiens mitochondrially encoded 12S ribosomal RNA (RNR1), ribosomal RNA
       https://www.ncbi.nlm.nih.gov/nuccore/NR_137294.1?report=fasta
       
       NR_137295.1 Homo sapiens mitochondrially encoded 16S ribosomal RNA (RNR2), ribosomal RNA
       https://www.ncbi.nlm.nih.gov/nuccore/NR_137295.1?report=fasta
    
   There is a copy fo UCSC tRNA and NCBI rRNA sequences to in *sequences* directory.
   
   Next Bowtie index will be built for the sequences, for example:
   
   ```bash
   # All contaminant sequeneces can be merged in a file, assume named contaimination.fa
   
   bowtie-build <contaimination.fa> <contaimination_bowtie_index_dir> 
   ```
   
   Then reads will be mapped to contaminant sequeneces, the unmapped reads will be kept for genome alignment, for example:
   
   ```bash
   bowtie -a --best --strata -S --seed 23 -p <threads> --chunkmbs 256 --norc --maqerr=60 --un <trimmed_unfiltered.fq> <contaimination_bowtie_index_dir> <(gzip -dc <trimmed.fq.gz>) <trimmed_filtered.sam>
   
   gzip <trimmed_unfiltered.fq>
   
   ```
   
#### 4. Aligning to reference genome

   We use START aligner. Firstly, we build STAR index for the genome.
   
   ```bash
   # sjdbOverhang is the longest read length - 1
   
   STAR --runThreadN <threads> --runMode genomeGenerate --genomeDir <star_ribo_dir> --genomeFastaFiles <ref_genome.fa> --sjdbGTFfile <annotation.gtf> --sjdbOverhang <sjdbOverhang>
   ```
   
   Then align to the reference genome, for example:
   
   ```bash
   SAM_ATTR="NH HI NM MD AS"
   ALIGNINTRON_MIN=20
   ALIGNINTRON_MAX=10000
   MISMATCH_MAX=1
   MISMATCH_NOVERL_MAX=0.04
   FILTER_TYPE=BySJout
   
   STAR --runThreadN <threads> --genomeDir <star_ribo_dir> \
        --readFilesIn <trimmed_unfiltered.fq.gz> --readFilesCommand zcat \
        --outReadsUnmapped Fastx --outFileNamePrefix <prefix> \
        --alignIntronMin <alignIntronMin> --alignIntronMax <alignIntronMax> --alignEndsType EndToEnd \
        --outFilterMismatchNmax <mismatch_max> --outFilterMismatchNoverLmax <mismatch_noverl_max>\
        --outFilterType <filter_type> --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMattributes <sam_attr> --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN <threads>
        
   # Only uniquely mapped reads are kept for ORF calling
   
   samtools view -b -q 255 <Aligned.sortedByCoord.out.bam> > <sample_q255.bam>
   samtools index <sample_q255.bam>
   ```

#### 5. P-site calling

   P-site offset is the distance from the 5’ or 3’ end of a ribosome-protected footprint (RPF) to the P-site of the ribosome that generated the footprint. Because the P-site is the site where peptidyl elongation occurs, read alignments from ribosome profiling are frequently mapped to their P-sites. We use plastid Python package to estimate it from Ribo-Seq data.

   Firstly, we create a protein coding gene annotation file (GTF format):
    
   ```bash
   # We keep protein coding leve 1 and 2, not seleno proteins
   
   cat <gencode_annotation.gtf> | awk '{if($18=="\"protein_coding\";" && $0~"level (1|2);" && $0!~"tag \"seleno\";" && $0!~"cds_end_NF" && $0!~"cds_start_NF"){print $0}}' | sort -k1,1 -k4,4n | bgzip > <protein_coding.gtf.gz>
   
   tabix -p gff <protein_coding.gtf.gz>
   ```

   Secondly, we run psite estimation and phasing estimation:
   
   ```bash
   metagene generate <metagene_analysis_dir> \
                     --landmark cds_start \
                     --annotation_files <protein_coding.gtf.gz>
   
   psite <metagene_analysis_dir>/protein_coding_rois.txt <psite_dir> \
         --min_length 25 \
         --max_length 35 \
         --require_upstream \
         --count_files <sample_q255.bam>

   phase_by_size <metagene_analysis_dir>/protein_coding_rois.txt <phasing_dir> \
                 --count_files <sample_q255.bam> \
                 --fiveprime_variable \
                 --offset <psite_dir>/p_offsets.txt \
                 --codon_buffer 5
   ```
   
   The output of for psite offset estimation is like: 
   
       length  p_offset
       25      7
       26      8
       27      9
       28      12
       29      12
       30      12
       31      12
       32      13
       33      13
       34      13
       35      13
       default 13
       
   The output of phasing is like:
   
       read_length     reads_counted   fraction_reads_counted  phase0  phase1  phase2
       25      19781   0.015942        0.332946        0.443304        0.223750
       26      26675   0.021498        0.331621        0.264105        0.404274
       27      47856   0.038568        0.588829        0.177094        0.234077
       28      97003   0.078176        0.461099        0.166923        0.371978
       29      216829  0.174745        0.607806        0.175069        0.217125 
       30      374353  0.301695        0.699292        0.086186        0.214522
       31      298826  0.240827        0.534723        0.052951        0.412327
       32      121240  0.097709        0.592610        0.329627        0.077763
       33      30141   0.024291        0.567732        0.313991        0.118277
       34      6519    0.005254        0.557908        0.312931        0.129161
       35      1609    0.001297        0.517091        0.324425        0.158484
       
   We can see that read length 29-31 have strong biased for phase0 or reading frame 1 (triplet periodicity, a feature of Ribo-Seq data)    

### RNA-Seq data processing

QC and adapter trimming are similar to Ribo-Seq data. RNA-seq data can be longer than 50bp (e.g. 75bp/100bp) and paired-end, when making STAR index, the sjdbOverhang needs to change accordingly (e.g. if the library is 100bp, the sjdbOverhang should be 100-1). One example to align paired-end RNA-seq using STAR:

```bash
nice -5 STAR --runThreadN <threads> --genomeDir <star_rna_100_pe_dir> \
                       --readFilesIn <R1_trimmed.fq.gz> <R2_trimmed.fq.gz> --readFilesCommand zcat \
                       --alignEndsType EndToEnd \
                       --outFilterMismatchNmax <mismatch_max> \
                       --outReadsUnmapped Fastx --outFileNamePrefix <prefix> \
                       --outSAMattributes <sam_attr> --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
                       --outBAMsortingThreadN <threads>
```

#### Transcript expression estimation

One additional step is to estimate transcript expression value (FPKM or TPM), only expressed transcripts (FPKM > 0.5) will be considered for ORF calling. We use StringTie to estimate FPKM values, for example:  

```bash
stringtie <sample_q255.bam>  -p <threads> -G <protein_coding.gtf> -eB -o expressed.gtf -A gene_abund.tab -C cov_refs.gtf
```

### ORF calling

The ORF calling workflow has the following steps:

#### 1. BAM filter

Given a range of read length (min_len, max_len) based on phasing estimation, the BAM file will be filtered, for example:

```bash
samtools view -h <ribo.bam> | awk -v minl="<min_len>" -v maxl="<max_len>" 'length($10) >= minl && length($10) <= maxl || $1 ~ /^@/' | samtools view -bS - > filtered.bam

samtools index filtered.bam
```

#### 2. Read count filter

Filter out ORFs with 0 read covered. bedtools is used to count the reads, for example:

```bash
parallel -j<threads> "bedtools coverage -counts -split -a <orfs_{}.smORFs.bed> -b <filtered.bam> > <orfs_{}.smORFs.bedtools.counts>" ::: ATG CTG TTG GTG

parallel "awk '\$13 > 0' <orfs_{}.smORFs.bedtools.counts> > <orfs_{}.smORFs.bedtools.greaterThanZeroReads.counts>" ::: ATG CTG TTG GTG
parallel "cut -f1-12 <orfs_{}.smORFs.bedtools.greaterThanZeroReads.counts> > <orfs_{}.smORFs.bedtools.greaterThanZeroReads.bed>" ::: ATG CTG TTG GTG
```

#### 3. RPF count filter

Not all reads are RPFs, plastid will extract all RPFs given p-site offset information, for example:

```bash
parallel "python script/get_count_vectors.py --annotation_files <orfs_{}.smORFs.bedtools.greaterThanZeroReads.bed> --annotation_format BED --count_files <filtered.bam> --min_length <min_len> --max_length <max_len> --fiveprime_variable --offset <psite_offset_file> --out_prefix {} ." ::: ATG CTG TTG GTG
```

#### 4. Transcript expression filter

We will only consider those smORFs whose host transcripts are expressed (by stringtie estimation). 

```bash
parallel "LC_ALL=C fgrep -f <stringytie_active_transcript_Ids.txt> <plastid_ORFs_{}_translated.txt> > <ORFs_{}_translated_expressed.txt>" ::: ATG CTG TTG GTG
```

#### 5. Label filter

According to their location on host transcript, smORFs will be classified as following:

Class | Description
----- | -----------
canonical | overlapping annotated CDS, have the same stop with annotated CDS
canonical_extended | overlapping annotated CDS, have the same stop with annotated CDS
canonical_truncated | overlapping annotated CDS, have the same stop with annotated CDS
five_prime | in upstream of annotated CDS
five_prime_overlap | in upstream of annotated CDS, overlapping annotated CDS
three_prime | in downstream of annotated CDS
three_prime_overlap | in downstream of annotated CDS, overlapping annotated CDS
within | in internal of annotated CDS, but in a different frame relative annotated CDS
noncoding | in non-coding genes or non-coding transcripts of coding genes

#### 6. ORFScore filter

ORFScore will be calculated and any smORFs with NA score and low coverage (< 0.1) will be discared, for example:

```bash
parallel "Rscript <ORFScore.R> <ORFs_{}_translated_expressed.txt> {}" ::: ATG CTG TTG GTG

cat <ORFScore_*TG.tsv> | sort -k2 -g | awk '{if($2!="NA") print}' > <ORFScore_all_noNA.tsv>
cat <ORFScore_all_noNA.tsv> | awk '{if($2 >0 && $6 >= 0.1) print}' > <ORFScore_all_filteredByCov.tsv>
```

#### 7. Region filter

If a smORF region is overlapping with annotated coding seqeunces (CDSs), estimate the propotion of signal (RPF_CDS/RPF_smORF) it obsorbs from the CDS, if it is greater than 1, the smORF is inside the CDSs and will be filtered out, for example:

```bash
Rscript ~/code/github/orf-discovery/script/RegionFilter.R . <cds.gtf> <filtered.bam> <threads>
```

#### 8. Ribosome release score (RRS) filter

Calculate RRS, filter out rna_ratio < 45 and RRS < 5 (cutoff based on [RRS paper](https://www.ncbi.nlm.nih.gov/pubmed/23810193))

```bash
Rscript script/RRS.R . <rna-seq_bam_dir>

tail -n +2 <RRS_out.tsv> | awk '$2 > 9 && $3 < 30 && $4 > 5' > <RRS_out_filtered.tsv>
```

#### 9. Nested filter

ORFs may have the same stop but different start, they are nested, among the nested ORFs, the one with ATG start codon has the highest priority, then the one with max ORFScore will be selected, for example:

```bash
Rscript script/NestedRegionFilter.R . <threads>
```

#### 10. FDR filter

ORFscore p-vals will be adjusted (q-val), and we set 0.01 as cutoff. In the meantime, we'd like keep the smORFs longer than 10 codons. A stringent constraint is that we expect to see RPFs covering codon 1 and 2, so the mean of RPF is greater than 0. For example:

```bash
Rscript script/ORFScore_padj.R <ORFScore_merged_PT_filtered.tsv> <all_withQval>

awk '($4+0) < 0.01' <ORFScore_all_withQval.tsv> > <ORFScore_fdr_filtered.tsv>

# Length filter >= 11aa including stop codon
cat <ORFScore_fdr_filtered.tsv> | cut -f1 | cut -f4,5 -d: | tr ":" "\t" | awk '{print ($2-$1+1)/3}' > <ORFScore_fdr_filtered_len.txt>
paste <ORFScore_fdr_filtered.tsv> <ORFScore_fdr_filtered_len.txt> | awk '$17 >= 11' | cut -f1-16 > <ORFScore_length_filtered.tsv>

# codon_1_2 mean should be greater than 0
awk '$10 > 0' ORFScore_length_filtered.tsv >ORFScore_codon_filtered.tsv
```

## Support

### email

Please report any issues or questions by creating a ticket, or by email to 
<fengyuan.hu@babraham.ac.uk>.
