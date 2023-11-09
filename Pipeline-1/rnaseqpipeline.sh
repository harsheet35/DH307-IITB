#!/bin/bash

# Step 1: Quality control using FastQC
fastqc ERR5021531_1.fastq.gz ERR5021531_2.fastq.gz -o fastqc_output

# Step 2: Trim adapters and low-quality bases using Trimmomatic
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 ERR5021531_1.fastq.gz ERR5021531_2.fastq.gz \
    output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
    ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

# Step 3: Alignment using HISAT2
hisat2 -q -x ~/Desktop/Harsheet_DH307/grch38/genome -1 output_R1_paired.fastq.gz -2 output_R2_paired.fastq.gz -S output_alignment.sam																																																																																																																																																																																																																																																																																																												

# Step 4: Convert SAM to BAM and sort
samtools view -bS output_alignment.sam | samtools sort -o output_alignment_sorted.bam

# Step 5: Index the sorted BAM file
samtools index output_alignment_sorted.bam

# Step 6: Assemble transcripts and quantify using StringTie
stringtie output_alignment_sorted.bam -G ~/Desktop/Harsheet_DH307/grch38/genome.gtf -o output_transcripts.gtf

# Step 7: Run differential expression analysis using DESeq2 in R
# Rscript deseq2.R

# End of RNA-Seq Pipeline
