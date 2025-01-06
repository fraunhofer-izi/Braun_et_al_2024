
## Workflow for the whole genome sequencing (WGS) data including variant calling, reverse engineering of the cilta-cel vector sequence, and STAR alignment for integration site analysis.

![fig_11_supp](https://github.com/user-attachments/assets/e426ef41-a77f-4f69-9a2d-84f63332aca9)

## STAR = Detection of candidates for integration loci of cilta-cel  
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

## viralintegration 
* viralintegration pipeline used (nf-core Pipeline 'viralintegration' https://nf-co.re/viralintegration/0.1.1/docs/output/)
* run screen for complete viral databank (no input for --viral_fasta necessary, default setting used)
* to checked for Epstein-Barr-Virus (EBV) integration only, specify viral_fasta in params.yaml with: 'viral_fasta ./data/EBV/GCF_002402265.1_ASM240226v1_genomic.fna'
	* The viral fasta reference for EBV (human gammaherpesvirus 4) was downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/.

Steps:
* 01_run_viralintegration.sh
