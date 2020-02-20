#!/bin/sh
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=25gb
#SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --job-name=align
#SBATCH --error=/work/soybean/mhapp95/AGRO932/AGRO932/log/align.error
#SBATCH --output=/work/soybean/mhapp95/AGRO932/AGRO932/log/align.out

cd /work/soybean/mhapp95/AGRO932/AGRO932/data

module load bwa samtools

###### use only if reference hasn't been indexed yet ######

# index the reference genome
#bwa index reference/ZeaMays_Chloroplast.fasta 

###########################################################


for i in {1..10}; do bwa mem reference/ZeaMays_Chloroplast.fasta raw/pop1_$i.read1.fq raw/pop1_$i.read2.fq | 
samtools view -bSh - > raw/pop1_$i.bam; 
done
# sort
for i in {1..10}; do samtools sort raw/pop1_$i.bam -o raw/pop1_sorted$i.bam; done
# index them
for i in {1..10}; do samtools index raw/pop1_sorted$i.bam; done

for i in {1..10}; do bwa mem reference/ZeaMays_Chloroplast.fasta raw/pop2_$i.read1.fq raw/pop2_$i.read2.fq | samtools view -bSh - > raw/pop2_$i.bam; 
done
# sort
for i in {1..10}; do samtools sort raw/pop2_$i.bam -o raw/pop2_sorted$i.bam; done
# index them
for i in {1..10}; do samtools index raw/pop2_sorted$i.bam; done


