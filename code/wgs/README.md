# STAR = Detection of candidates for integration loci of cilta-cel  
* use parameters from arriba to find chimeric reads
* enable mutlti mappers in order to not exlude LTR regions
* do not filter secondary reads in IGV
* filter for alignment score >= 100
* filter for unique mappers: NH:i:1 and HI:i:1
* requirements for integration site candidates:
	* true hg38-Carvykti chimeric read
	* differing positions of similar candidates have to be within a range of 10 bases to be recognized as the same candidate
	* assignment flanking reads have to be within 200 bases of the integration site to be counted as support for the candidate
	* minimum of 2 reads have to support a candidate to nominate a candidate
* requirements scripted in R

Steps:
* download GRCh38.primary_assembly.genome.fa and GENCODE v46 annotation -> ./data/
* generate custom reference with cilta-cel -> ./data/
* Perform steps 01_ till 07_ 

# viralintegration_EBV = Investigation of Epstein-Barr-Virus (EBV) integration 
* viralintegration pipeline used (nf-core Pipeline 'viralintegration' https://nf-co.re/viralintegration/0.1.1/docs/output/) for detection of EBV to exclude EBV as the cause for tumorigenesis

Steps:
* 01_run_viralintegration.sh