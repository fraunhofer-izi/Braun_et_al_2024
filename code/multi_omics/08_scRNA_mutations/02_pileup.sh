#!/bin/bash

module load BCFtools

POS="chr4:105235572-105235572,chr1:64860244-64860244,chr11:36400995-36400995,chr22:28687938-28687938"

# Find all *.bai files in current directory and pileup at relevant genomic positions
find ./data/scRNA_BAMs_sorted/split_bams -type f -name "*.bam" | while read -r bam_file; do
  bcftools mpileup --fasta-ref ./data/GRCh38.primary_assembly.genome.fa --annotate FORMAT/AD -r "$POS" "$bam_file" > "${bam_file%.bam}_mpileup.txt"
done


