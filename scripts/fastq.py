##Last updated: 7/8/26

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that downloads fastq files and joins them into one file via FASTQ-DUMP on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('srr', metavar='SRAs', nargs='+', help='SRR names to download for individual strain')
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-s", metavar="STRAIN", dest="strain", help="Strain name for naming files", required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)

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
#SBATCH --job-name=fastq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=""")
fastq.write(straindir)
fastq.write("""/%x_%j.out
""")
if args.email is not None:
	fastq.write("""#SBATCH --mail-user=""")
	fastq.write(args.email)
	fastq.write("""
#SBATCH --mail-type=FAIL""")
fastq.write('\n')
fastq.write('\n')

##write load modules
fastq.write("""module purge
ml SRA-Toolkit/3.0.10-gompi-2023a
ml SAMtools/1.18-GCC-12.3.0""")
fastq.write('\n')
fastq.write('\n')

##change directories
fastq.write("".join(['cd ',userdir,'/', strain, '\n', '\n']))

##write fastq-dump command
for srr in SRRlist:
	fastq.write("".join(['fastq-dump ', '--split-3 ', '--gzip ', srr, '\n']))

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
