##Last updated:	10/2/25

##import system
import sys
import subprocess
import os
import argparse
import numpy as np
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import datetime
import string
from decimal import Decimal
parser = argparse.ArgumentParser(description='Script for removing duplicate CNAGs on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
new_seqs=f"{userdir}/Orthogroup_Sequences/"
os.makedirs(os.path.dirname(new_seqs),exist_ok=True)

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
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

##Go through all the files
for file in os.listdir(orthodir):
	filename=("".join([orthodir, '/', file]))
	new_filename=("".join([new_seqs, file]))
	orthogroup = os.path.splitext(file)[0]
	CNAGs=[]
	all_seqs=[]

	for record in SeqIO.parse(filename,"fasta"):
		id=record.id
		if "CNAG" in id:
			no_prot=id.split("-")[0]
			if no_prot not in CNAGs:
				CNAGs.append(no_prot)
		else:
			all_seqs.append(record)

	for CNAG in CNAGs:
		longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
		for record in SeqIO.parse(filename,"fasta"):
			if CNAG in record.id:
				if len(record)>len(longest):
					longest=record
		if longest:
			rename=SeqRecord(longest.seq, longest.id,"","")
			all_seqs.append(rename)

	sorted_records=sorted(all_seqs,key=lambda r:r.id)
	print(all_seqs)

	SeqIO.write(sorted_records,new_filename,"fasta")
