##Last updated: 7/7/26

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
parser = argparse.ArgumentParser(description='Script to complete miniprot.py data',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-m", metavar="MINI_DIR", dest="minidir", help="miniprot directory created by miniprot.py", required=True)

args = parser.parse_args()

##Check files exist
if os.path.isdir(args.minidir) and os.access(args.minidir, os.R_OK):
	minidir=(os.path.abspath(args.minidir))
else:
	sys.exit("ERROR: miniprot directory is missing or not readable!!")

##file definitions
userdir = os.getcwd()
outdir=f"{userdir}/mini_final/"
error_out=f"{outdir}short_read_missing.csv"

##Go through each directory
missing_orthos={}
genome_count={}
mapped_count={}
file_count=0
if not os.path.isdir(outdir):
	os.makedirs(os.path.dirname(outdir),exist_ok=False)

for genome in os.listdir(minidir):
	genome_dir=f"{minidir}/{genome}"
	header=genome[0:10].replace("_","-")
	missing_hits=f"{genome_dir}/{genome}_missing_orthos.csv"
	current_orthos=f"{genome_dir}/{header}.csv"
	mapped_ortho_file=f"{genome_dir}/{genome}_mapped_orthos.csv"
	called_orthos=[]
	mapped_orthos=[]
	if os.path.getsize(missing_hits)==0:
		print("".join(["long_mini.py did not run correctly for ",genome]))
	if os.path.isfile(missing_hits) and os.path.isfile(current_orthos) and os.path.isfile(mapped_ortho_file):
		file_count+=1
		with open(current_orthos,"r") as fin:
			for line in fin:
				ortho=line.split(",")[0]
				if "M" in ortho:
					ortho=(line.split(",")[0]).split("M")[0]
				if ortho not in called_orthos:
					called_orthos.append(ortho)
					if ortho in genome_count:
						genome_count[ortho].append(genome)
					else:
						genome_count[ortho]=[genome]
		with open(missing_hits,"r") as fin2:
			for line in fin2:
				ortholog=line.split("\n")[0]
				if ortholog not in called_orthos:
					if ortholog in missing_orthos:
						missing_orthos[ortholog].append(genome)
					else:
						missing_orthos[ortholog]=[genome]
		with open(mapped_ortho_file,"r") as fin3:
			for line in fin3:
				ortho=line.split("\n")[0]
				if ortho in called_orthos:
					if ortho not in mapped_orthos:
						mapped_orthos.append(ortho)
						if ortho in mapped_count:
							mapped_count[ortho].append(genome)
						else:
							mapped_count[ortho]=[genome]

with open(error_out,"w") as f:
	f.write("".join(['missing_ortholog,false_positive_rate,false_negative_rate,sensitivity,specificity,accuracy\n']))
	for ortholog in genome_count:
		if ortholog not in mapped_count:
			mapped_count[ortholog]=[]
		if ortholog not in missing_orthos:
			missing_orthos[ortholog]=[]
		true_positive=len(set(mapped_count[ortholog]).intersection(set(genome_count[ortholog])))
		true_negative=file_count-len(list(set(genome_count[ortholog])-set(missing_orthos[ortholog])))
		false_positive=len(list(set(genome_count[ortholog])-set(mapped_count[ortholog])))
		false_negative=len(missing_orthos[ortholog])
		accuracy=(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)
		try:
			false_positive_rate=false_positive/(false_positive+true_negative)
		except ZeroDivisionError:
			false_positive_rate="NA"
		try:
			false_negative_rate=false_negative/(false_negative+true_positive)
		except ZeroDivisionError:
			false_negative_rate="NA"
		try:
			sensitivity=true_positive/(true_positive+false_positive)
		except ZeroDivisionError:
			sensitivity="NA"
		try:
			specificity=true_negative/(true_negative+false_positive)
		except ZeroDivisionError:
			specificity="NA"
		#print("".join([ortholog,",",str(true_positive),",",str(true_negative),",",
		#		str(false_positive),",",str(false_negative)]))
		f.write("".join([ortholog,",",str(false_positive_rate),",",str(false_negative_rate),",",
				str(sensitivity),",",str(specificity),",",str(accuracy),"\n"]))

