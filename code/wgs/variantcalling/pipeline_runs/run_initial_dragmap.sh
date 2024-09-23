#!/bin/bash

# this is entirely optional, it just prevents downloading the singularity containers multiple times
export NXF_SINGULARITY_CACHEDIR="/some/accessible/cachedir/"

module load Singularity/3.8.7
module load Nextflow/23.10.0

mkdir -p initial_dragmap/
cd initial_dragmap/

# Run Sarek pipeline on slurm/singularity

nextflow run nf-core/sarek -r 3.4.3 \
  --input ../samplesheet_plus_reseq.csv \
  -c ../custom.config \
  -profile slurm,singularity,rib20limited \
  -resume \
  --max_cpus 32 \
  --max_memory '240.GB' \
  --validationFailUnrecognisedParams \
  --outdir sarek_GATK_hg38_plus_Carvykti_dragmap/ \
  --aligner dragmap \
  --skip_tools baserecalibrator \
  --dragmap false \
  --save_reference \
  --genome GATK.GRCh38 \
  --fasta /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.fasta \
  --fasta_fai /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.fasta.fai \
  --three_prime_clip_r1 1 \
  --three_prime_clip_r2 1 \
  --igenomes_base /path/to/igenomes/ \
  --vep_cache /path/to/vep_cache \
  --vep_include_fasta \
  --dict /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.dict \
  --intervals /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed \
  --tools deepvariant,manta,strelka,mutect2,vep

# disabled because of incompatibility with dragmap for some arcane reason, see https://github.com/Illumina/DRAGMAP/issues/11
#  --trim_fastq \
#  --trim_nextseq 10 \
