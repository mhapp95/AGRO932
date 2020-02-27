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
angsd -bam raw/bam.txt -doSaf 1 -anc reference/ZeaMays_Chloroplast.fasta -GL 1 -out raw/out 
realSFS raw/out.saf.idx > raw/out.sfs

###### calculate theta #####
angsd -bam raw/bam.txt -out raw/out -doThetas 1 -doSaf 1 -pest raw/out.sfs -anc reference/ZeaMays_Chloroplast.fasta -GL 1
thetaStat print raw/out.thetas.idx > raw/theta.txt
thetaStat do_stat raw/out.thetas.idx -win 5000 -step 1000  -outnames thetasWindow

##### calculate Fst #####

#ls raw/pop1_sorted*.bam > raw/pop1bam.txt
#ls raw/pop2_sorted*.bam > raw/pop2bam.txt

#angsd -b raw/pop1bam.txt -anc reference/ZeaMays_Chloroplast.fasta -out raw/pop1 -dosaf 1 -gl 1
#angsd -b raw/pop2bam.txt -anc reference/ZeaMays_Chloroplast.fasta -out raw/pop2 -dosaf 1 -gl 1

# calculate the 2dsfs prior
#realSFS raw/pop1.saf.idx raw/pop2.saf.idx > raw/pop1.pop2.ml

# prepare the fst for easy window analysis etc
realSFS fst index raw/pop1.saf.idx raw/pop2.saf.idx -sfs raw/pop1.pop2.ml -fstout raw/out

# get the global estimate
realSFS fst stats raw/out.fst.idx > raw/fst_global.txt

# below is not tested that much, but seems to work
realSFS fst stats2 raw/out.fst.idx -win 500 -step 100 > raw/fst_win.txt
