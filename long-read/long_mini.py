##Last updated: 6/16/26

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
parser = argparse.ArgumentParser(description='Code that runs miniprot on long read assemblies to get missing sequences on GACRC Slurm',
epilog='Written by Claire Reichardt')

parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY_DIR", dest="asbdir", help="Directory containing list of long read assembilies", required=True)
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
required.add_argument("-p",dest='protdir', metavar='PROT_DIR', help='Directory of protein fasta files',required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
LOGFILE = f"{args.orthodir}/Log.txt"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
	logfile = (os.path.abspath(LOGFILE))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(args.asbdir) and os.access(args.asbdir, os.R_OK):
	asbdir=(os.path.abspath(args.asbdir))
else:
	sys.exit("ERROR: Assembly directory is missing or not readable!!")

if os.path.isdir(args.protdir) and os.access(args.protdir, os.R_OK):
	protdir=(os.path.abspath(args.protdir))
else:
	sys.exit("ERROR: protein directory is missing or not readable!!")

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
				if "cryptococcus_" in infile:
					genomename=infile.split("cryptococcus_")[1]
				else:
					genomename=infile
				if genomename not in genomes:
					no_n=genomename.split("\n")[0]
					genomes.append(no_n)
					no_fa=os.path.splitext(no_n)[0]
					header=no_fa[0:10].replace('_','-')
					headers.append(header)
if len(genomes)>0:
	print("".join(['Found ', str(len(genomes)), ' genomes in : ', orthodir]))
else:
	sys.exit("ERROR: improperly formatted logfile in orthofinder directory")


##Write miniprot files for each sequence
long_genomes=[]
os.makedirs(os.path.dirname(minidir),exist_ok=True)
for file in os.listdir(asbdir):
	genome_name=os.path.splitext(file)[0]
	header=genome_name[0:10].replace('_','-')
	protein_file=f"{protdir}/{header}.fa"
	if not os.path.isfile(protein_file):
		sys.exit("".join(['ERROR: Could not find protein file for ', genome_name,' in ', protdir]))
	mini_gdir=f"{minidir}{genome_name}/"
	shellfile=f"{mini_gdir}minishell.sh"
	assembly_file=f"{asbdir}/{file}"
	index_file=f"{mini_gdir}{genome_name}.mpi"
	gff_file=f"{mini_gdir}{genome_name}_aln.gff"
	missing_file=f"{mini_gdir}{genome_name}_missing_orthos.fa"
	missing_gff=f"{mini_gdir}{genome_name}_missing.gff"
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
	if args.email is not None:
		shell.write("""#SBATCH --mail-user=""")
		shell.write(args.email)
		shell.write("""
#SBATCH --mail-type=FAIL""")
	shell.write('\n')
	shell.write('\n')

	#write load modules
	shell.write("""module purge
ml miniprot/0.13-GCCcore-12.3.0
ml Python/3.11.3-GCCcore-12.3.0
ml Biopython/1.83-foss-2023a
ml GenomeTools/1.6.5-GCC-12.3.0""")
	shell.write("".join(['\n', '\n']))

	##change directories
	shell.write("".join(['cd ', mini_gdir, '\n\n']))

	#index the genome
	shell.write("".join(['miniprot -t 8 -d ',index_file, ' ',assembly_file,'\n']))

	#get gff for each search
	shell.write("".join(['miniprot -t 8 -I --gff --gff-only --outn 1 -u --gff-delim M ',index_file,' ',
	protein_file,' > ', gff_file,'\n\n']))

	#execute get-orthos.py script
	shell.write("".join(['python ',userdir,'/get-orthos.py ',genome_name,' ', orthodir,'\n\n']))

	#get gff for missing orthologs
	shell.write("".join(['miniprot -t 8 -I --gff --gff-only --outc 0.8 --gff-delim M --trans ',
	index_file,' ', missing_file,' > ', missing_gff,'\n\n']))

	#execute check-proteins.py script
	shell.write("".join(['python ',userdir,'/mini-overlap.py ',genome_name,'\n\n']))

	#help with submission
	shell.close()
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)

