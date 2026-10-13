##Last updated: 7/24/24

##import system
import sys
import subprocess
import os
import argparse
import fnmatch
parser = argparse.ArgumentParser(description='Code that runs AAFTF for genome assemblies on GACRC Slurm',
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
shellfiles=[]
userdir = os.getcwd()
jobsfile=f"{userdir}/assemblies/jobs.txt"
os.makedirs(os.path.dirname(jobsfile),exist_ok=True)

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
		shellfile = f"{userdir}/assemblies/{strain}/{strain}_assembly.sh"
		os.makedirs(os.path.dirname(shellfile), exist_ok=True)
		shellfiles.append(shellfile)
		straindir = f"{userdir}/assemblies/{strain}"
		forward_trim = f"{strain}_1P.fastq.gz"
		reverse_trim = f"{strain}_2P.fastq.gz"
		mito_out = f"{strain}.mito.fasta"
		forward_filter = f"{strain}_filtered_1.fastq.gz"
		reverse_filter = f"{strain}_filtered_2.fastq.gz"
		spades_out = f"{strain}.spades.fasta"
		vecscreen_out = f"{strain}.vecscreen.fasta"
		sourpurge_out = f"{strain}.sourpurge.fasta"
		rmdup_out = f"{strain}.rmdup.fasta"
		polish_out = f"{strain}.polish.fasta"
		sort_out = f"{strain}.sort.fasta"
		assess_out = f"{strain}.assess.txt"

		##create a file for each header
		shell = open(shellfile, "w")
		shell.write("""#!/bin/bash
#SBATCH --job-name=aaftf
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=48:00:00
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
		#Trim
		shell.write("".join(['AAFTF trim -l ', forward, ' -r ', reverse, '\n']))
		#Mito
		shell.write("".join(['AAFTF mito -l ', forward_trim, ' -r ', reverse_trim, ' -o ', mito_out, '\n']))
		#Filter
		shell.write("".join(['AAFTF filter -l ', forward_trim, ' -r ', reverse_trim, ' -o ', strain, '\n']))
		#Assemble
		shell.write("".join(['AAFTF assemble -l ', forward_filter, ' -r ', reverse_filter, ' -c 24 -o ', spades_out, '\n']))
		#Vecscreen
			##SKIP THIS STEP BECAUSE IT IS TAKING OUT BIG CHUNKS OF THE GENOME
			##Think this is related to mitochondrial genes that have relocated to the genome^
		#shell.write("".join(['AAFTF vecscreen -i ', spades_out, ' -c 16 -o ', vecscreen_out, '\n']))

		#Sourpurge
		shell.write("".join(['AAFTF sourpurge -i ', spades_out, ' -o ', sourpurge_out, ' -c 24 --phylum ', args.phylum, '\n']))
		#Remove duplicates
		shell.write("".join(['AAFTF rmdup -i ', sourpurge_out, ' -o ', rmdup_out, '\n']))
		#Polish
		shell.write("".join(['AAFTF polish -i ', rmdup_out, ' -l ', forward, ' -r ', reverse, ' -o ', polish_out, '\n']))
		#Sort
		shell.write("".join(['AAFTF sort -i ', polish_out, ' -o ', sort_out, '\n']))
		#Assess
		shell.write("".join(['AAFTF assess -i ', sort_out, ' -r ', assess_out, '\n']))
		shell.close()

#write execute line in jobsfile
sorted_shells=sorted(shellfiles)
with open(jobsfile, "w") as f:
	for shellfile in sorted_shells:
		f.write("".join(['sbatch ', shellfile, '\n']))

##print to execute shell file
print("".join(['Found ', str(count), ' paired sequence files for assembly']))
if count ==0:
        sys.exit("".join(['ERROR: No valid sequence files found in: ', srrdir,
        ' Please check if input files are formatted correctly and end in: .fastq.gz']))
else:
	pass

##Only submit first 100:
with open(jobsfile, "r") as fin:
	lines=fin.readlines()
batchdir = f"{userdir}/assemblies/sbatchexec/"
os.makedirs(os.path.dirname(batchdir), exist_ok=True)
countlines = 0
countbatch = 0
start = 0
end = len(lines)
step = 90

#Count jobs:
for i in range(start,end):
	x = i
	countlines +=1
print("".join(['Found ', str(countlines), ' jobs in: ', jobsfile, ' for submission to slurm.']))

#Make large batch files
for i in range(start,end,step):
	x=i
	countbatch +=1
	batchname = ("".join([batchdir, 'batch_', str(countbatch),'.sh']))
	with open (batchname, 'w') as f:
		f.write("""#!/bin/bash
#SBATCH --job-name=aaftf
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH --dependency=singleton
#SBATCH --output=""")
		f.write("".join([batchdir, 'batch_', str(countbatch),'.out']))
		f.write('\n')
		f.write('\n')
		lotsolines=("".join(lines[x:x+step]))
		f.write(lotsolines)
	f.close()

	#Help with submission
	if args.submit is True:
		execute=("".join(['sbatch ', batchname, '\n']))
		os.system(execute)
	else:
		print("".join(['--submit flag was not used. For submission of jobs to slurm please use: sbatch ',
		batchname]))
