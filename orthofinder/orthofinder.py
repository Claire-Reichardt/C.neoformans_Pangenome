##Last updated:	8/4/25

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that makes a shell to run Orthofinder GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
parser.add_argument('proteindir', metavar='PROTEIN_DIRECTORY', help='Protein files (.protein.fa) folder from Funannotate output')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
diamonddir = f"{userdir}/Diamond/"
os.makedirs(os.path.dirname(diamonddir), exist_ok=True)
shellfile = f"{userdir}/oshell.sh"
os.makedirs(os.path.dirname(shellfile), exist_ok=True)

##Check files exist
PROTEINDIR = args.proteindir
if os.path.isdir(PROTEINDIR) and os.access(PROTEINDIR, os.R_OK):
	proteindir = (os.path.abspath(args.proteindir))
else:
	sys.exit("ERROR: Protein sequence directory is missing or not readable!!")

##check dependent scripts are in the same folder
diamondscript = f"{userdir}/diamond.py"
if os.access(diamondscript, os.R_OK):
	pass
else:
	sys.exit(''.join(["ERROR: Dependency script diamond.py not found in current working directory:", userdir,
	". Please move diamond.py to current directory to successfully run all steps of OrthoFinder."]))
finishscript =f"{userdir}/finish.py"
if os.access(finishscript, os.R_OK):
	pass
else:
        sys.exit(''.join(["ERROR: Dependency script finish.py not found in current working directory:", userdir,
        ". Please move diamond.py to current directory to successfully run all steps of OrthoFinder."]))

##print files
count = 0
for fp in os.listdir(proteindir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count +=1
	elif ext==".fasta":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any protein files in: ', proteindir, '. Please check provided files are correctly formatted.']))
else:
	print("".join(['Found ', str(count), ' protein files in: ', proteindir]))

##create a file for each header
shell = open(shellfile, "w")
shell.write("""#!/bin/bash
#SBATCH --job-name=ortho
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=""")
shell.write(userdir)
shell.write("""/%x.out
""")
if args.email is not None:
	shell.write("""#SBATCH --mail-user=""")
	shell.write(args.email)
	shell.write("""
#SBATCH --mail-type=FAIL""")
shell.write('\n')
shell.write('\n')

##write load modules
shell.write("""module purge
ml OrthoFinder/2.5.5-foss-2023a
ml Python/3.11.3-GCCcore-12.3.0""")
shell.write('\n')
shell.write('\n')

##Run orthofinder
ortho = ['orthofinder -f ', proteindir, ' -op -o ', userdir, '/Orthofinder_diamond', '\n']
printortho = ''.join(ortho)
shell.write(printortho)

##get diamond execution scripts
diamond = ["grep ^diamond ortho.out | grep -v 'commands that must be run' > ", diamonddir, "jobs.txt", '\n']
printdiamond = ''.join(diamond)
shell.write(printdiamond)

##execute diamond submodule
shell.write("python diamond.py")

##print to execute shell file
shell.close()
if args.submit is True:
	execute = ("".join(['sbatch ', shellfile, '\n']))
	os.system(execute)
