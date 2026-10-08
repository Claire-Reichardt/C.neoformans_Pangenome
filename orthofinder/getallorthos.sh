#!/bin/bash
#SBATCH --job-name=getallorthos
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/cer92177/orthofinder/%x_%j.out
#SBATCH --mail-user=claire.reichardt@uga.edu
#SBATCH --mail-type=END,FAIL

ml Python/3.11.3-GCCcore-12.3.0
ml Biopython/1.83-foss-2023a

python getallorthos.py -o /scratch/cer92177/orthofinder/Orthofinder_diamond/Results_Aug05/
