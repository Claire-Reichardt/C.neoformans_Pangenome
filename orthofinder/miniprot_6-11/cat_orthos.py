##Last updated: 2/12/26

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
parser = argparse.ArgumentParser(description='Script to combine orthologs based on mapping data',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
required.add_argument('-mo', metavar='ORTHOLOG_OVERLAP',dest="ortholog_overlap", help='.csv file containing overlapping information to combine orthologs')

args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
LOGFILE = f"{args.orthodir}/Log.txt"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
	logfile = (os.path.abspath(LOGFILE))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isfile(args.ortholog_overlap) and os.access(args.ortholog_overlap,os.R_OK):
	overlap=os.path.abspath(args.ortholog_overlap)
else:
	sys.exit("ERROR: Overlap file is missing or not readable!!")

##file definitions
userdir = os.getcwd()
combodir=f"{userdir}/combine/"

##Go through and get orthologs that overlap
if not os.path.isdir(combodir):
	os.makedirs(os.path.dirname(combodir),exist_ok=False)

	#Copy sequence directory to edit outside of original files
	for file in os.listdir(orthodir):
		filename=f"{orthodir}/{file}"
		shutil.copy2(filename,combodir)

	#Create a dictionary of orthologs to combine and go from smallest to largest ortholog
	combine_dict={}
	deleted_files={}
	with open(overlap,'r') as fin:
		for line in fin:
			combine_decision=line.split(",")[7]
			if combine_decision=="TRUE":
				ortho1=line.split(",")[0]
				ortho2=line.split(",")[1]
				if ortho2 in deleted_files:
					upper_ortho=deleted_files[ortho2]
					if upper_ortho in combine_dict:
						combine_dict[upper_ortho].append(ortho2)
						combine_dict[upper_ortho].append(ortho1)
					else:
						combine_dict[upper_ortho]=[ortho2]
						combine_dict[upper_ortho].append(ortho1)
					print("".join(['Combining ', ortho2,' with ',upper_ortho, ' instead of ',
						ortho1,' because of prior combinations.']))
					deleted_files[ortho2]=upper_ortho
					deleted_files[ortho1]=upper_ortho

				if ortho2 not in deleted_files:
					deleted_files[ortho2]=ortho1
					if ortho1 in combine_dict:
						combine_dict[ortho1].append(ortho2)
					else:
						combine_dict[ortho1]=[ortho2]

	sorted_combo_dict=dict(sorted(combine_dict.items(),reverse=True))

	#Go through dictionary, remove files and combine orthologs
	for ortho1,ortho2s in sorted_combo_dict.items():
		unique_ortho2s=[]
		for ortho2 in ortho2s:
			if ortho2 not in unique_ortho2s:
				unique_ortho2s.append(ortho2)
				new_seqs=[]
				ortho1_file=f"{combodir}{ortho1}.fa"
				ortho2_file=f"{combodir}{ortho2}.fa"
				if os.path.isfile(ortho1_file):
					if os.path.isfile(ortho2_file):
						for record in SeqIO.parse(ortho1_file,"fasta"):
							new_seqs.append(record)
						for record2 in SeqIO.parse(ortho2_file,"fasta"):
							new_seqs.append(record2)
						sorted_seqs=sorted(new_seqs,key=lambda r:r.id)
						SeqIO.write(sorted_seqs,ortho1_file,"fasta")

						os.remove(ortho2_file)
						print("".join(['Combined ', ortho1,' and ', ortho2]))

					else:
						print("".join(['cannot find ',ortho2]))

				else:
					print("".join([ortho1," has already been combined, please resolve this error."]))


else:
	print("".join([combodir, ' already exists. Please move directory to run this script.']))
