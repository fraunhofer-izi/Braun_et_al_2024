#!/bin/bash

ml STAR

STAR --runMode genomeGenerate --runThreadN 30 --genomeDir STAR_index_Carvykti/ --genomeFastaFiles ./data/hg38_Carvykti_vector.fasta --sjdbGTFfile ./data/gencode.v46.primary_assembly.basic.annotation.gtf sjdbOverhang 149 > LOG_ascii 2> ERR_ascii
