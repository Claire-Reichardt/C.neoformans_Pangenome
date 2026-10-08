##Last updated: 2/18/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
parser = argparse.ArgumentParser(description='submodule of miniprot.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
parser.add_argument('ortho_dir')
args = parser.parse_args()

##file definitions
dir = os.getcwd()

##get all missing orthologs
directory_header=args.dir[0:10].replace("_","-")
current_orthos=f"{dir}/{directory_header}.csv"
if not os.path.isfile(current_orthos):
	current_orthologs={}
	for file in os.listdir(args.ortho_dir):
		orthogroup=os.path.splitext(file)[0]
		filename=f"{args.ortho_dir}/{file}"
		for record in SeqIO.parse(filename,"fasta"):
			genome = (record.id).split("_")[0]
			if directory_header == genome:
				if record.description not in current_orthologs:
					current_orthologs[record.description]=orthogroup
	#write to csv
	with open(current_orthos,'w') as outfile:
		for prot,ortho in current_orthologs.items():
			outfile.write("".join([ortho,',',prot,'\n']))

	if len(current_orthologs)==0:
		print("".join(['Unable to find any orthologs for header: ', directory_header,
				'Please check this is the correct protein heading for genome: ',
					args.dir]))
else:
	print("".join(["Found current orthologs under ",current_orthos," Skipping get-orthos.py"]))
