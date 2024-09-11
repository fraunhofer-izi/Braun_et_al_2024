#!/bin/bash
#SBATCH --job-name=PB_CARpos_CD5dim
#SBATCH --output=%x_%j.txt
#SBATCH --error=%x_%j_gpu_error.txt
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time 72:00:00
#SBATCH --cpus-per-task 30
#SBATCH --mem 240000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user nhu-nguyen.pham@izi.fraunhofer.de

ml STAR

OLDFASTQ="./data/*fastq" #first-run S2
NEWFASTQ="./data/*fastq" #re-sequenced S2
INDEX="STAR_index_Carvykti/"

STAR --runMode alignReads \
	--runThreadN 30 \
	--genomeDir ${INDEX} \
	--readFilesCommand zcat \
	--readFilesIn ${OLDFASTQ}PB_CARpos_CD5dim_S2_R1_001.fastq.gz,${NEWFASTQ}P2251_CARpos_CD5dim_S1_R1_001.fastq.gz ${OLDFASTQ}PB_CARpos_CD5dim_S2_R2_001.fastq.gz,${NEWFASTQ}P2251_CARpos_CD5dim_S1_R2_001.fastq.gz \
	--outSAMattrRGline ID:run1 , ID:run2 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outReadsUnmapped Fastx \
	--twopassMode Basic \
	--limitSjdbInsertNsj 5000000 \
	--outFilterMultimapNmax 50 \
	--peOverlapNbasesMin 10 \
	--alignSplicedMateMapLminOverLmate 0.5 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 6 \
	--chimOutType Junctions WithinBAM HardClip \
	--chimJunctionOverhangMin 10 \
	--chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 \
	--chimScoreSeparation 1 \
	--chimSegmentReadGapMax 3 \
	--chimMultimapNmax 50 \
	--winAnchorMultimapNmax 100
