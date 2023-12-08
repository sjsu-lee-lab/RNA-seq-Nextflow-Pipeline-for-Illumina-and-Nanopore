#!/bin/bash 

#SBATCH --job-name=nextflow-test 
#SBATCH --output=%x.o%j 
#SBATCH --error=%x.e%j  

cd ~/RNASeq/
./nextflow run rnaseq.nf -ansi-log false -resume $1