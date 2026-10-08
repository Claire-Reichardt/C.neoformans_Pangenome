##Last updated:	5/28/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description='Code to extract CNAG annotation from protein.fa on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('protfile', metavar='PROTEIN_FILE', help='H99 protein file to extract annotations from')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outfile = f"{userdir}/cnags.csv"

##Check files exist
PROTFILE = args.protfile
if os.path.isfile(PROTFILE) and os.access(PROTFILE, os.R_OK):
	protfile = (os.path.abspath(args.protfile))
else:
	sys.exit("ERROR: Input protein file is missing or not readable!")

#Search for CNAG
for record in SeqIO.parse(protfile, "fasta"):
	if "CNAG" in record.id:
		with open(outfile, "a") as f:
			f.write(record.id)
			splitter = record.description.split("|")
			for item in splitter:
				if "gene_product" in item:
					f.write(",")
					product_only=item.split("=")[1]
					f.write(product_only)
					f.write("\n")
