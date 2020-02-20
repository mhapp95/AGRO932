#!/bin/sh
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=25gb
#SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --job-name=angsd
#SBATCH --error=/work/soybean/mhapp95/AGRO932/AGRO932/log/angsd.error
#SBATCH --output=/work/soybean/mhapp95/AGRO932/AGRO932/log/angsd.out

cd /work/soybean/mhapp95/AGRO932/AGRO932/data/

##### write bam to txt file #####
ls raw/*sorted*.bam > raw/bam.txt
angsd -bam raw/bam.txt -doSaf 1 -anc reference/GlyMax_CompleteMitochondria.fasta -GL 1  -out raw/out 
realSFS raw/out.saf.idx > raw/out.sfs

###### calculate theta #####
angsd -bam raw/bam.txt -out raw/out -doThetas 1 -doSaf 1 -pest raw/out.sfs -anc reference/GlyMax_CompleteMitochondria.fasta -GL 1
thetaStat print raw/out.thetas.idx > raw/theta.txt
