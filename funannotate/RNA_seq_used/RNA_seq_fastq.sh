#!/bin/bash
#SBATCH --job-name=fastq                        # Job name
#SBATCH --partition=batch                       # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                              # 1 task (process) for below commands
#SBATCH --cpus-per-task=8                      # CPU core count per task, #by default 1 CPU core per task
#SBATCH --mem=20G                              # Memory per node (4GB); by default using M as unit
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/cer92177/fun_6-9/2ofeach/%x_%j.out                      # Standard output log, e.g., testBowtie2_12345.out
#SBATCH --mail-user=Claire.Reichardt@uga.edu    # Where to send mail
#SBATCH --mail-type=FAIL                        # Mail events (BEGIN, END, FAIL, ALL)

ml SRA-Toolkit/3.0.1-centos_linux64
ml SAMtools/1.16.1-GCC-11.3.0

cd /scratch/cer92177/fun_6-9/2ofeach

fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip SRR554166
fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip SRR554139
fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip SRR554154
fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip SRR554185
fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --gzip SRR554182
cat *_1.fastq.gz > 2ofeach_1.fq.gz
cat *_2.fastq.gz > 2ofeach_2.fq.gz
rm *fastq.gz
