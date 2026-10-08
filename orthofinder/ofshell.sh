#!/bin/bash
#SBATCH --job-name=ortho
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=3-00:00:00
#SBATCH --dependency=singleton
#SBATCH --output=/scratch/cer92177/orthofinder/%x_%j.out


module purge
ml OrthoFinder/2.5.5-foss-2023a
ml Python/3.11.3-GCCcore-12.3.0
ulimit -Hn
ulimit -Sn

orthofinder --fewer-files -b /scratch/cer92177/orthofinder/Orthofinder_diamond/Results_Aug04
