#!/bin/bash
#SBATCH --job-name=ortho
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/cer92177/orthofinder/%x.out
#SBATCH --mail-user=Claire.Reichardt@uga.edu
#SBATCH --mail-type=FAIL

module purge
ml OrthoFinder/2.5.5-foss-2023a
ml Python/3.11.3-GCCcore-12.3.0

orthofinder -f /scratch/cer92177/orthofinder/split_protein_files_7-23 -op -o /scratch/cer92177/orthofinder/Orthofinder_diamond
grep ^diamond ortho.out | grep -v 'commands that must be run' > /scratch/cer92177/orthofinder/Diamond/jobs.txt
python diamond.py