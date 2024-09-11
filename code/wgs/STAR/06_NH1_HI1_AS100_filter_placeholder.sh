#!/bin/bash

ml SAMtools

SAMPLE=$(pwd | cut -d '/' -f 9) #this is to give the samples an appropiate namen and can be changed

samtools view -e '[AS]>=100' Aligned.sortedByCoord.out.bam |\
       grep Carvykti_vector |\
       awk -F " " '$12 == "NH:i:1"' |\
       awk -F " " '$13 == "HI:i:1"' >\
       ${SAMPLE}.Carvykti_vector.NH1_HI1_AS100_filtered.sam

grep Carvykti_vector Chimeric.out.junction > ${SAMPLE}.Chimeric.Carvykti_vector.txt

cut -f 10 ${SAMPLE}.Chimeric.Carvykti_vector.txt > ${SAMPLE}.Chimeric.Carvykti_vector.readnames.txt

grep -f ${SAMPLE}.Chimeric.Carvykti_vector.readnames.txt ${SAMPLE}.Carvykti_vector.NH1_HI1_AS100_filtered.sam |\
       sort |\
       uniq |\
       cut -f 1 >\
       ${SAMPLE}.filtered.reads.txt

grep -f ${SAMPLE}.filtered.reads.txt ${SAMPLE}.Chimeric.Carvykti_vector.txt > ${SAMPLE}.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt

awk -F " " '$7 == 0' ${SAMPLE}.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt > ${SAMPLE}.Chimeric.Carvykti_vector.filtered.Type0.txt
awk -F " " '$7 == -1' ${SAMPLE}.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt > ${SAMPLE}.Chimeric.Carvykti_vector.filtered.flanking.txt

cut -f 10 ${SAMPLE}.Chimeric.Carvykti_vector.filtered.Type0.txt > ${SAMPLE}.Chimeric.Carvykti_vector.filtered.Type0_readnames.txt
cut -f 10 ${SAMPLE}.Chimeric.Carvykti_vector.filtered.flanking.txt > ${SAMPLE}.Chimeric.Carvykti_vector.filtered.flanking_readnames.txt
