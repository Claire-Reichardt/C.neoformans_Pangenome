##Last updated: 3/15/24

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that downloads fastq files and joins them into one file via FASTQ-DUMP on GACRC Slurm - for FUNANNOTATE formatting!!',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('srr', metavar='SRAs', nargs='+', help='SRR names to download for individual strain')
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-s", metavar="STRAIN", dest="strain", help="Strain name for naming files", required=True)
required.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm progress", required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
strain = args.strain
SRRlist = args.srr
straindir = f"{userdir}/{strain}"
fastqfile = f"{userdir}/{strain}/{strain}_fastq.sh"
os.makedirs(os.path.dirname(fastqfile), exist_ok=True)

##create a header for each file
fastq = open(fastqfile, 'w')
fastq.write ("""#!/bin/bash
#SBATCH --job-name=fastq                        # Job name
#SBATCH --partition=batch                       # Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1                              # 1 task (process) for below commands
#SBATCH --cpus-per-task=24                      # CPU core count per task, #by default 1 CPU core per task
#SBATCH --mem=100G                              # Memory per node (4GB); by default using M as unit
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=""")
fastq.write(straindir)
fastq.write("""/%x_%j.out                      # Standard output log, e.g., testBowtie2_12345.out
""")
fastq.write("""#SBATCH --mail-user=""")
fastq.write(args.email)
fastq.write("""    # Where to send mail
#SBATCH --mail-type=FAIL                        # Mail events (BEGIN, END, FAIL, ALL)"""
)
fastq.write('\n')
fastq.write('\n')

##write load modules
fastq.write("""ml SRA-Toolkit/3.0.1-centos_linux64
ml SAMtools/1.16.1-GCC-11.3.0""")
fastq.write('\n')
fastq.write('\n')

##change directories
fastq.write("".join(['cd ', straindir, '\n', '\n']))

##write fastq-dump command
for srr in SRRlist:
	fastq.write("".join(
	['fastq-dump ', "--split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' ", '--gzip ', srr, '\n']))

##join forward and reverse files
fastq.write("".join(['cat *_1.fastq.gz > ', strain, '_1.fq.gz', '\n']))
fastq.write("".join(['cat *_2.fastq.gz > ', strain, '_2.fq.gz', '\n']))

##remove raw files
fastq.write("".join(['rm *fastq.gz', '\n']))

##print to execute shell file
fastq.close()
if args.submit is True:
	execute = ("".join(['sbatch ', fastqfile, '\n']))
	os.system(execute)
