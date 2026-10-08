##Last updated:	8/5/25

##import system
import sys
import subprocess
import os
from datetime import datetime
import argparse
parser = argparse.ArgumentParser(description='Submodule for running Orthofinder GACRC Slurm',
epilog='Written by Claire Reichardt')

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outlog = f"{userdir}/ortho.out"
orthodir = f"{userdir}/Orthofinder_diamond/"
diamonddir = f"{userdir}/Diamond/"

##get date run
with open (outlog, "r") as f:
	for line in f:
		if "Starting OrthoFinder" in line:
			date=(line[0:10])
formatdate = datetime.strptime(date,'%Y-%m-%d')
monthandday = formatdate.strftime("%b%d")
resultsdate = f"{orthodir}Results_{monthandday}"
#check this is a folder
if os.path.isdir(resultsdate) and os.access(resultsdate, os.R_OK):
	print(''.join(["Found previous Diamond results in: ", resultsdate]))
else:
	sys.exit("ERROR: Could not OrthoFinder results folder. Please check previous steps have run correctly")

##Check all Diamond searches have completed successfully
for fp in os.listdir(diamonddir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".sh":
		diamondrun = os.path.splitext(fp)[0]
		outfile = ''.join([diamonddir, diamondrun, ".out"])
		if os.path.isfile(outfile) and os.access(outfile, os.R_OK):
			pass
		else:
			print(outfile)
			sys.exit(''.join(["Missing Diamond runs for: ", diamondrun,
			". Please check previous steps have run correctly."]))


##Run last step of Orthofinder
shellfile = f"{userdir}/ofshell.sh"
os.makedirs(os.path.dirname(shellfile), exist_ok=True)

#create a file for each header
shell = open(shellfile, "w")
shell.write("""#!/bin/bash
#SBATCH --job-name=ortho
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=3-00:00:00
#SBATCH --dependency=singleton
#SBATCH --output=""")
shell.write(userdir)
shell.write("""/%x_%j.out
""")
shell.write('\n')
shell.write('\n')

#write load modules
shell.write("""module purge
ml OrthoFinder/2.5.5-foss-2023a
ml Python/3.11.3-GCCcore-12.3.0
ulimit -Hn
ulimit -Sn""")
shell.write('\n')
shell.write('\n')

#Run orthofinder
ortho = ['orthofinder --fewer-files -b ', resultsdate,'\n']
printortho = ''.join(ortho)
shell.write(printortho)

##print to execute shell file
shell.close()
execute = ("".join(['sbatch ', shellfile, '\n']))
os.system(execute)
