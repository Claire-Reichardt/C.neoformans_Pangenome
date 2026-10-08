##Last updated:	5/30/24

##import system
import sys
import subprocess
import os
import argparse
import fnmatch
parser = argparse.ArgumentParser(description='Code that reruns AAFTF for genome assemblies after spAdes step on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-s", metavar="SRRDIRECTORY", dest="srrdir", help="Directory containing both forward and reverse reads of all sequences formatted in .fq.gz", required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required.add_argument("-p", metavar="PHYLUM", dest="phylum", help="Phylum of fungi species for sourmash to avoid parsing (ex. Basidiomycota)", required=True)
args = parser.parse_args()

##Check files exist
SRR = args.srrdir
if os.path.isdir(SRR) and os.access(SRR, os.R_OK):
	srrdir = (os.path.abspath(args.srrdir))
else:
	sys.exit("ERROR: SRR directory is missing or not readable!!")

##Write shell files for every strain
count = 0
for fp in os.listdir(srrdir):
	if fnmatch.fnmatch(fp, '*_1.fastq.gz'):
		count +=1
		strain = fp.replace('_1.fastq.gz','')
		forward = f"{srrdir}/{strain}_1.fastq.gz"
		if os.path.isfile(forward) and os.access(forward, os.R_OK):
			forward = os.path.abspath(forward)
		reverse = f"{srrdir}/{strain}_2.fastq.gz"
		if os.path.isfile(reverse) and os.access(reverse, os.R_OK):
			reverse= os.path.abspath(reverse)
		else:
			sys.exit(''.join(["ERROR:Could not find matching reverse reads for:",
			forward, ". Please make sure forward(_1) and reverse(_2) files are named the same."]))

		##file definitions
		userdir = os.getcwd()
		shellfile = f"{userdir}/assemblies/{strain}/{strain}_assembly.sh"
		os.makedirs(os.path.dirname(shellfile), exist_ok=True)
		jobsfile = f"{userdir}/assemblies/jobs.txt"
		os.makedirs(os.path.dirname(jobsfile), exist_ok=True)
		straindir = f"{userdir}/assemblies/{strain}"
		spades_out = f"{srrdir}/{strain}.spades.fasta"
		sourpurge_out = f"{strain}.sourpurge.fasta"
		rmdup_out = f"{strain}.rmdup.fasta"
		polish_out = f"{strain}.polish.fasta"
		sort_out = f"{strain}.sort.fasta"
		assess_out = f"{strain}.assess.txt"

		##create a file for each header
		shell = open(shellfile, "w")
		shell.write("""#!/bin/bash
#SBATCH --job-name=aaftf
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=""")
		shell.write(straindir)
		shell.write("""/%x_%j.out
""")
		if args.email is not None:
			shell.write("""#SBATCH --mail-user=""")
			shell.write(args.email)
			shell.write('\n')
			shell.write("""#SBATCH --mail-type=FAIL""")
		shell.write('\n')
		shell.write('\n')

		##write load modules
		shell.write("""module purge
ml SRA-Toolkit/3.0.1-centos_linux64
ml SAMtools/1.16.1-GCC-11.3.0
export AAFTF_DB=~/lib/AAFTF_DB
source /home/${USER}/.bashrc
source activate aaftf""")
		shell.write('\n')
		shell.write('\n')

		##change directories
		shell.write("".join(['cd ', straindir, '\n']))
		shell.write('\n')

		##write AAFTF functions
		#Sourpurge
		shell.write("".join(['AAFTF sourpurge -i ', spades_out, ' -o ', sourpurge_out, ' -c 12 --phylum ', args.phylum, '\n']))

		#Remove duplicates
		shell.write("".join(['AAFTF rmdup -i ', sourpurge_out, ' -o ', rmdup_out, ' -c 12', '\n']))

		#Polish
		shell.write("".join(['AAFTF polish -i ', rmdup_out, ' -l ', forward, ' -r ', reverse,
		' -o ', polish_out, ' -c 12', '\n']))

		#Sort
		shell.write("".join(['AAFTF sort -i ', polish_out, ' -o ', sort_out, '\n']))

		#Assess
		shell.write("".join(['AAFTF assess -i ', sort_out, ' -r ', assess_out,'\n']))
		shell.close()

		#Help with submission
		if args.submit is True:
			execute=("".join(['sbatch ', shellfile, '\n']))
			os.system(execute)
		else:
			print("".join(['--submit flag was not used. For submission of jobs to slurm please use: sbatch ',
			shellfile]))

##print to execute shell file
print("".join(['Found ', str(count), ' paired sequence files for assembly']))
if count ==0:
        sys.exit("".join(['ERROR: No valid sequence files found in: ', srrdir,
        ' Please check if input files are formatted correctly and end in: .fq.gz']))
else:
	pass
