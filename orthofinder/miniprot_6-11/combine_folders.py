##Last updated: 2/17/26

##import system
import sys
import subprocess
import os
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
from pathlib import Path
import shutil
parser = argparse.ArgumentParser(description='Script to combine orthogroup data',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-m", metavar="MINI_DIR", dest="minidir", help="miniprot directory created by miniprot.py", required=True)
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)

args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
LOGFILE = f"{args.orthodir}/Log.txt"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
	logfile = (os.path.abspath(LOGFILE))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

##file definitions
userdir = os.getcwd()
minidir=os.path.abspath(args.minidir)
complete_dir=f"{userdir}/final_orthogroups/"

##Go through each
if not os.path.isdir(complete_dir):
	os.makedirs(os.path.dirname(complete_dir),exist_ok=False)
	for file in os.listdir(orthodir):
		filename=f"{orthodir}/{file}"
		minifile=f"{minidir}/{file}"
		if os.path.isfile(minifile):
			new_file=f"{complete_dir}/{file}"
			new_seqs=[]
			for record in SeqIO.parse(filename,"fasta"):
				new_seqs.append(record)
			for record2 in SeqIO.parse(minifile,"fasta"):
				new_seqs.append(record2)
			sorted_seqs=sorted(new_seqs,key=lambda r:r.id)
			SeqIO.write(sorted_seqs,new_file,"fasta")
		else:
			shutil.copy2(filename,complete_dir)
