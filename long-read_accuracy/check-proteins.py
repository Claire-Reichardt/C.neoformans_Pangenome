##Last updated: 6/1/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
parser = argparse.ArgumentParser(description='submodule of long_mini.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
args = parser.parse_args()

##file definitions
dir = os.getcwd()

##get all missing orthologs
directory_header=args.dir[0:10].replace("_","-")
current_orthos=f"{dir}/{directory_header}.csv"
gff_file=f"{dir}/{args.dir}_aln.gff"
if os.path.isfile(current_orthos) and os.path.isfile(gff_file):
	annotated_proteins=[]
	found_proteins=[]
	with open(current_orthos) as fin:
		for line in fin:
			if "OG" in line:
				protein=((line.split(",")[1]).split("\n")[0]).split(" ")[0]
			else:
				protein=(line.split(",")[1]).split("\n")[0]
			annotated_proteins.append(protein)
	with open(gff_file) as fin:
		for line in fin:
			if "##" not in line:
				feature=line.split("\t")[2]
				if feature=="mRNA":
					att=(line.split("\t")[8]).split("\n")[0]
					protein=((att.split(";")[0]).split("ID=")[1]).split("M")[0]
					ident=float((att.split(";")[2]).split("Identity=")[1])
					if ident>0.85:
						found_proteins.append(protein)
	missing_prots=list(set(annotated_proteins) - set(found_proteins))
	print(f"Proteins not found: {missing_prots}")
	print(str(len(missing_prots)))

