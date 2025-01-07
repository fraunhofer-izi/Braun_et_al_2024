#!/bin/bash
#SBATCH --job-name=Carvykti_merged_WGS_EBV
#SBATCH --output=%x_%j.txt
#SBATCH --error=%x_%j_gpu_error.txt
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time 120:00:00
#SBATCH --cpus-per-task 30
#SBATCH --mem 240000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user nhu-nguyen.do@izi.fraunhofer.de

ml Singularity
ml Nextflow/23.10.0

nextflow run nf-core/viralintegration \
       	-profile singularity \
	-params-file params.yaml \
	-c custom.config \
	-r 0.1.1 -resume