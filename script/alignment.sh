
bowtie -a --best --strata -m 100 --seed 23 -p 8 --chunkmbs 256 ../../../../ref/ensembl/mouse/transcriptome/Mus_musculus.GRCm38.cdna.protein_coding_transcripts -S <(gzip -dc ncrna_unmapped.fasta) | samtools view -F 0xC -bS - >transcriptome_mapped_seed23.bam; samtools sort transcriptome_mapped_seed23.bam transcriptome_mapped_seed23_sorted; samtools index transcriptome_mapped_seed23_sorted.bam

bowtie2 -k 100 --seed 23 -p 8 -x ../../../../ref/ensembl/mouse/transcriptome/Mus_musculus.GRCm38.cdna.protein_coding_transcripts -U <(gzip -dc ncrna_unmapped.fasta) | samtools view -F 0xC -bS - >transcriptome_mapped_bowtie2.bam

