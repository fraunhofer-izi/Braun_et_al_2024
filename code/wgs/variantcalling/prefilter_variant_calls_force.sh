#!/bin/bash

module load BCFtools/1.14-GCC-11.2.0

SCRIPTDIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPTDIR

# location to write prefiltered output to
OUTDIR="variant_calls_prefiltered/final"
mkdir -p $OUTDIR
OUTFILE="${OUTDIR}/Koeln_1.mutect2.filtered_VEP.normalized.vcf.gz"

# locations of sarek pipeline output, relative to this script
FORCE_CALL_RESULTS="pipeline_runs/forcecall_dragmap/sarek_GATK_hg38_plus_Carvykti_dragmap_force/annotation/mutect2/Koeln_1/Koeln_1.mutect2.filtered_VEP.ann.vcf.gz"

# use bcftools norm to normalize indels and split multiallelic sites into multiple rows as for the initial results for comparison

echo "Normalizing joint call file: Koeln_1.mutect2.filtered_VEP.ann.vcf.gz"

bcftools norm --threads 8 \
  --fasta-ref "/mnt/ribolution/user_worktmp/florian.peter.grosse/scratch/igenomes/Homo_sapiens_plus_Carvykti/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38_plus_Carvykti.fasta" \
  -m - \
  -o $OUTFILE $FORCE_CALL_RESULTS
tabix $OUTFILE
