#!/bin/bash

# index vcf if necessary
cp ../merged_forced_alleles.vcf.gz
tabix -f merged_forced_alleles.vcf.gz

# this is entirely optional, it just prevents downloading the singularity containers multiple times
export NXF_SINGULARITY_CACHEDIR="/some/accessible/cachedir/"

module load Singularity/3.8.7
module load Nextflow/23.10.0

mkdir -p forcecall_dragmap/
cd forcecall_dragmap/

# Run Sarek pipeline on slurm/singularity

nextflow run nf-core/sarek -r 3.4.3 \
  -c ../forcecall_custom.config \
  -profile slurm,singularity,rib20limited \
  -resume \
  --max_cpus 32 \
  --max_memory '240.GB' \
  --input /path/to/pipeline_runs/initial_dragmap/sarek_GATK_hg38_plus_Carvykti_dragmap/csv/markduplicates_no_table.csv \
  --step variant_calling \
  --joint_mutect2 \
  --validationFailUnrecognisedParams \
  --outdir sarek_GATK_hg38_plus_Carvykti_dragmap_force/ \
  --genome GATK.GRCh38 \
  --fasta /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.fasta \
  --fasta_fai /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.fasta.fai \
  --igenomes_base /path/to/igenomes/ \
  --vep_cache /path/to/vep_cache \
  --vep_include_fasta \
  --dict /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.dict \
  --intervals /path/to/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed \
  --tools mutect2,vep,ascat,msisensorpro
