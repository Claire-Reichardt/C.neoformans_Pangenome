##Last updated: 8/25/26

##import system
import sys
import subprocess
import os
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
parser = argparse.ArgumentParser(description='Code that runs miniprot to get see accessory hits on GACRC Slurm',
epilog='Written by Claire Reichardt')

parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY_DIR", dest="asbdir", help="Directory containing list of assembilies used for gene prediction", required=True)
required.add_argument('-d',dest='accdir', metavar='ACCESSORY_DIRECTORY', help='Directory of accessory ortholog fasta files)', required=True)
args = parser.parse_args()


##Check files exist
ORTHODIR = f"{args.accdir}"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Accessory orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(args.asbdir) and os.access(args.asbdir, os.R_OK):
	asbdir=(os.path.abspath(args.asbdir))
else:
	sys.exit("ERROR: Assembly directory is missing or not readable!!")

##file definitions
userdir = os.getcwd()
minidir=f"{userdir}/miniprot/"

##print files
count = 0
for fp in os.listdir(orthodir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any Orthogroup files in: ', orthodir, '. Please check OrthoFinder ran completely in provided directory.']))
else:
	print("".join(['Found ', str(count), ' Orthogroup files in: ', orthodir]))

##Now write miniprot files for each sequence
os.makedirs(os.path.dirname(minidir),exist_ok=True)
for file in os.listdir(asbdir):
	genome_name=os.path.splitext(file)[0]
	mini_gdir=f"{minidir}{genome_name}/"
	shellfile=f"{mini_gdir}minishell.sh"
	assembly_file=f"{asbdir}/{file}"
	index_file=f"{mini_gdir}{genome_name}.mpi"

	os.makedirs(os.path.dirname(shellfile),exist_ok=True)
	shell=open(shellfile,"w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=miniprot
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --output=""")
	shell.write(f"{mini_gdir}")
	shell.write("""%x_%j.out
""")
	shell.write('\n')
	shell.write('\n')

	#write load modules
	shell.write("""module purge
ml miniprot/0.13-GCCcore-12.3.0""")
	shell.write("".join(['\n', '\n']))

	##change directories
	shell.write("".join(['cd ', mini_gdir, '\n\n']))

	#index the genome
	shell.write("".join(['miniprot -t 8 -d ',index_file, ' ',assembly_file,'\n']))

	for orthogroup in os.listdir(orthodir):
		ortho=os.path.splitext(orthogroup)[0]
		gff_file=f"{mini_gdir}{genome_name}_{ortho}.gff3"
		prot_dir=f"{orthodir}/{orthogroup}"

		#get gff and protein for each search
		shell.write("".join(['miniprot -t 8 -I --outc 0.5 --gff --gff-only --gff-delim M --trans ',index_file,
		' ',prot_dir,' > ', gff_file,'\n']))

	shell.close()
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)

