##Last updated: 6/9/25

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that runs ragtag on assemblies from aaftf on GACRC Slurm',
epilog = 'Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm progress", default=None)
required.add_argument("-a", metavar="ASSEMBLIES", dest="assemblies", help="Directory with fasta assembly sequences", required=True)
required.add_argument("-r", metavar="REFERENCE", dest="ref", help="Reference genome for scaffolding contigs to", required=True)


args = parser.parse_args()

##file definitions
userdir = os.getcwd()
ragtag_dir=f"{userdir}/ragtag_assemblies/"
os.makedirs(os.path.dirname(ragtag_dir), exist_ok=True)

##check assembly file is valid
ASSEMBLYPATH = args.assemblies
if os.path.isdir(ASSEMBLYPATH) and os.access(ASSEMBLYPATH, os.R_OK):
    assemblypath= (os.path.abspath(args.assemblies))
else:
    sys.exit("ERROR: Assembly directory is missing or not readable!")

REFERENCE=args.ref
if os.path.isfile(REFERENCE) and os.access(REFERENCE, os.R_OK):
	reffile=(os.path.abspath(args.ref))
else:
	sys.exit("ERROR: Reference file is missing or not readable!")

#Create a shell for each protein
for file in os.listdir(assemblypath):
	filename=os.path.join(assemblypath,file)
	if ".fa" in file:
		strain=file.split(".fa")[0]
	straindir = f"{userdir}/{strain}"
	shellfile = f"{userdir}/{strain}/{strain}_rshell.sh"
	os.makedirs(os.path.dirname(shellfile), exist_ok=True)

	shell = open(shellfile, 'w')
	shell.write ("""#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=""")
	shell.write(straindir)
	shell.write("""/%x_%j.out
""")
	shell.write("""#SBATCH --mail-user=""")
	if args.email is not None:
		shell.write(args.email)
		shell.write("\n")
		shell.write("""#SBATCH --mail-type=FAIL""")
	shell.write('\n')
	shell.write('\n')

	##write load modules
	shell.write("""module purge
ml RagTag/2.1.0""")
	shell.write('\n')
	shell.write('\n')

	##change directories
	shell.write("".join(['cd ', straindir, '\n', '\n']))

	##write BUSCO functions
	shell.write("".join(['ragtag.py scaffold ', reffile,' ',filename,
	' -o ', straindir, "/",strain, '_ragtag -r -t 8 -u',
	'\n']))

	##copy final file to directory
	shell.write("".join(['mv ',straindir,'/',strain,'_ragtag/ragtag.scaffold.fasta ',
	ragtag_dir,strain, '.fa','\n']))

	##print to execute shell file
	shell.close()
	if args.submit is True:
		execute = ("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)
