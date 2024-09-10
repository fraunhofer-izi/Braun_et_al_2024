#!/bin/bash

module load SAMtools

# Find all *.bam files in current directory and all subdirectories and index
find ./data/scRNA_BAMs_sorted/split_bams -type f -name "*.bam" | while read -r bam_file; do
  samtools index "$bam_file"
done
