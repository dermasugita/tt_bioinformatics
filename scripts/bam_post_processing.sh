#!/bin/bash
samtools sort ../bam/Normal_chr17_part.bam -o ../bam/Normal_chr17_part.sorted.bam
samtools sort ../bam/Tumor_chr17_part.bam -o ../bam/Tumor_chr17_part.sorted.bam
samtools index ../bam/Normal_chr17_part.sorted.bam
samtools index ../bam/Tumor_chr17_part.sorted.bam
