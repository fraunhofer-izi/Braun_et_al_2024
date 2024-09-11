#!/bin/bash

ml STAR

STAR --runMode alignReads --runThreadN 30 --genomeDir STAR_index_Carvykti/ --readFilesCommand zcat --readFilesIn ./data/Apherese_CAR-_CD5pos_S1_R1_001.fastq.gz ./data/Apherese_CAR-_CD5pos_S1_R2_001.fastq.gz --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outReadsUnmapped Fastx --twopassMode Basic --limitSjdbInsertNsj 5000000 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 6 --chimOutType Junctions WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --winAnchorMultimapNmax 100
