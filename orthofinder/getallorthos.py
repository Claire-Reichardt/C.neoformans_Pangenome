##Last updated:	2/13/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser(description='Code to extract strain-specific protein sequences from Ortholog groups on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o', metavar='ORTHOLOG_DIRECTORY', dest = 'orthodir',required=True,
help='Orthologous Sequence directory (ex.Orthofinder_diamond/Results_Apr22)  from OrthoFinder')

args = parser.parse_args()

##file definitions
userdir = os.getcwd()

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

##Get headers to search from from Orthologues dir
with open(logfile, 'r') as fin:
	lines=fin.readlines()
search = []

for line in lines:
	#extract lines of protein files used
	if not "#" in line:
		if ":" in line:
			infile=line.split(": ")[1]
			ext=os.path.splitext(infile)[-1].lower()
			if ".fa" in ext or ".fasta" in ext:
				file=infile
				if "_" in file:
					genomename= file.split("_")[1]
				else:
					genomename=file
				#genomebase = os.path.splitext(genomename)[0]
				name = os.path.splitext(genomename)[0]
				header = name[0:10].replace('_', '-')
				search.append(header)
				outfile = f"{args.orthodir}/all_orthos/{name}.csv"
				if not os.path.isfile(outfile):
					os.makedirs(os.path.dirname(outfile), exist_ok=True)

					#Search for genome headers
					for file in os.listdir(orthodir):
						filename=("".join([orthodir, '/', file]))
						orthogroup = os.path.splitext(file)[0]
						for record in SeqIO.parse(filename, "fasta"):
							id=record.id
							genome = id.split("_")[0]
							if header == genome:
								with open(outfile, "a") as f:
									f.write(orthogroup)
									f.write(",")
									f.write(record.id)
									f.write("\n")


