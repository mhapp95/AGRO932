#!/bin/sh
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=65gb
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --job-name=wgsim
#SBATCH --error=/work/soybean/mhapp95/AGRO932/AGRO932/log/wgsim.error
#SBATCH --output=/work/soybean/mhapp95/AGRO932/AGRO932/log/wgsim.out

cd /work/soybean/mhapp95/AGRO932/AGRO932/data

for i in {1..10}
do
   wgsim reference/ZeaMays_Chloroplast.fasta -N 3500 -1 150 -2 150 -r 0.03 -e 0.01 -R 0 -X 0 raw/pop1_$i.read1.fq raw/pop1_$i.read2.fq
done

for i in {1..10}
do
   wgsim reference/ZeaMays_Chloroplast.fasta -N 3500 -1 150 -2 150 -r 0.03 -e 0.01 -R 0 -X 0 raw/pop2_$i.read1.fq raw/pop2_$i.read2.fq
done
