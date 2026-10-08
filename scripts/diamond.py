##Last updated:	6/7/24

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Submodule for Orthofinder execution on GACRC Slurm',
epilog='Written by Claire Reichardt')

##file definitions
userdir = os.getcwd()
diamonddir = f"{userdir}/Diamond/"
os.makedirs(os.path.dirname(diamonddir), exist_ok=True)
shellfile = f"{userdir}/oshell.sh"
os.makedirs(os.path.dirname(shellfile), exist_ok=True)
jobsfile=f"{diamonddir}jobs.txt"

##make shells
with open(jobsfile,'r') as fin:
	lines = fin.readlines()
count = 0
start = 0
end = len(lines)
step = 100
for i in range(start, end, step):
	x = i
	count +=1
	filename = ("".join([diamonddir, 'diamond_', str(count), '.sh']))
	with open(filename, 'w') as f:
		f.write("""#!/bin/bash
#SBATCH --job-name=ortho
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=""")
		f.write("".join([diamonddir, 'diamond_', str(count),'.out']))
		f.write('\n')
		f.write('\n')
		f.write("""module purge
ml OrthoFinder/2.5.4-foss-2022a
ml Python/3.10.4-GCCcore-11.3.0""")
		f.write('\n')
		f.write('\n')
		lotsolines=("".join(lines[x:x+step]))
		f.write(lotsolines)
	f.close()

##change last file to execute finish.py script
lastfile=("".join([diamonddir, 'diamond_', str(count), '.sh']))
if os.path.isfile(lastfile) and os.access(lastfile, os.R_OK):
	with open(lastfile, "a") as f:
		f.write('\n')
		f.write("python finish.py")
	f.close()

##submit all scripts to bash
diamondfiles=[]
for filename in os.listdir(diamonddir):
	ext=os.path.splitext(filename)[-1]
	if ext==".sh":
		base=os.path.splitext(filename)[0]
		number_only=base.split("_")[-1]
		diamondfiles.append(int(number_only))

sorted_files=sorted(diamondfiles, key=int)
for filename in sorted_files:
	execute = ("".join(['sbatch ', diamonddir,'diamond_', str(filename),'.sh', '\n']))
	os.system(execute)
