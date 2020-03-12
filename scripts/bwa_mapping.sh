#!/bin/bash
# Now you may start mapping paired-end reads to the reference sequence!
# Modify -t option according to your environment
bwa mem -M -t 2 ../reference/chr17.fa ../fastq/Normal_chr17_part_R1.fastq.gz ../fastq/Normal_chr17_part_R2.fastq.gz | samtools view -bS - > ../bam/Normal_chr17_part.bam
bwa mem -M -t 2 ../reference/chr17.fa ../fastq/Tumor_chr17_part_R1.fastq.gz ../fastq/Tumor_chr17_part_R2.fastq.gz | samtools view -bS - > ../bam/Tumor_chr17_part.bam
