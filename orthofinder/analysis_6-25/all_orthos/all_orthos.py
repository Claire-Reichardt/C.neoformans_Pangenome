##Last updated: 6/25/26

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
parser = argparse.ArgumentParser(description='Code that runs collects protein information for all genomes GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
LOGFILE = f"{args.orthodir}/Log.txt"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
	logfile = (os.path.abspath(LOGFILE))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

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

genomes=[]
headers=[]
with open(logfile, 'r') as fin:
	lines=fin.readlines()
for line in lines:
	if not "#" in line:
		if ":" in line:
			infile=line.split(": ")[1]
			ext=os.path.splitext(infile)[-1].lower()
			if ".fa" in ext or ".fasta" in ext:
				if "_" in infile:
					genomename=infile.split("_")[1]
				else:
					genomename=infile
				if genomename not in genomes:
					no_n=genomename.split("\n")[0]
					no_fa=os.path.splitext(no_n)[0]
					genomes.append(no_fa)
					header=no_fa[0:10].replace('_','-')
					headers.append(header)
if len(genomes)>0:
	print("".join(['Found ', str(len(genomes)), ' genomes in : ', orthodir]))
else:
	sys.exit("ERROR: improperly formatted logfile in orthofinder directory")

##file definitions
userdir = os.getcwd()
ortho_path=os.path.abspath(args.orthodir)
alldir=f"{ortho_path}/all_orthos/"
shell_dir=f"{alldir}shell_files/"

##Now write miniprot files for each sequence
os.makedirs(os.path.dirname(alldir),exist_ok=True)
for genome in genomes:
	genome_name=genome
	all_gdir=f"{alldir}{genome_name}/"
	shellfile=f"{shell_dir}{genome_name}_all-shell.sh"
	os.makedirs(os.path.dirname(shellfile),exist_ok=True)
	shell=open(shellfile,"w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=get_orthos
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=""")
	shell.write(f"{shell_dir}{genome_name}")
	shell.write("""_%j.out
""")
	if args.email is not None:
		shell.write("""#SBATCH --mail-user=""")
		shell.write(args.email)
		shell.write("""
#SBATCH --mail-type=FAIL""")
	shell.write('\n')
	shell.write('\n')

	#write load modules
	shell.write("""module purge
ml Python/3.11.3-GCCcore-12.3.0
ml Biopython/1.83-foss-2023a""")
	shell.write("".join(['\n', '\n']))

	#execute get-orthos.py script
	shell.write("".join(['python ',userdir,'/get-orthos.py ',genome_name,' ', orthodir,'\n\n']))


	shell.close()
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)
