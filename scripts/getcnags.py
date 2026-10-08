##Last updated:	4/25/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description='Code to extract CNAG from Ortholog groups on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHOLOG_DIRECTORY', help='Orthologous Sequence directory (ex.Orthofinder_diamond/Results_Apr22_1/Orthogroup_Sequences)  from OrthoFinder')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outfile = f"{userdir}/fasearch.tsv"

##Check files exist
ORTHODIR = args.orthodir
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(args.orthodir))
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

#Search for CNAG
for file in os.listdir(orthodir):
	filename=("".join([orthodir, '/', file]))
	orthogroup = os.path.splitext(file)[0]
	for record in SeqIO.parse(filename, "fasta"):
		if "CNAG" in record.id:
			with open(outfile, "a") as f:
				f.write(orthogroup)
				f.write("\t")
				f.write(record.id)
				f.write("\n")
f.close()
