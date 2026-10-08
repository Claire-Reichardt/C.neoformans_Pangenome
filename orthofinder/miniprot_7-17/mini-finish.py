##Last updated: 4/9/26

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
outdir=f"{userdir}/mini_final"
gff_dir=f"{outdir}/gff_files/"
new_orthos=f"{outdir}/new_orthogroups/"
overlap_out=f"{outdir}/overlapping_orthologs.csv"

##Go thorugh each directory
overlap_orthos={}
if not os.path.isdir(gff_dir):
	os.makedirs(os.path.dirname(gff_dir),exist_ok=False)
	if not os.path.isdir(new_orthos):
		os.makedirs(os.path.dirname(new_orthos),exist_ok=False)

		for genome in os.listdir(minidir):
			prot_count=0
			genome_dir=f"{minidir}/{genome}"
			header=genome[0:10].replace("_","-")
			current_orthos=f"{genome_dir}/{header}.csv"
			final_gff=f"{genome_dir}/{genome}.gff3"
			new_prots=f"{genome_dir}/{genome}_newprots.fa"
			overlap_file=f"{genome_dir}/{genome}_overlapping-orthologs.csv"

			if os.path.getsize(current_orthos) !=0:
				if not os.path.isfile(f"{gff_dir}{genome}.gff3"):
					shutil.copy2(final_gff,gff_dir)
				if os.path.isfile(new_prots):
					for record in SeqIO.parse(new_prots,"fasta"):
						orthogroup=((record.description).split(" ")[1]).split("M")[0]
						ortho_file=f"{new_orthos}{orthogroup}.fa"
						with open(ortho_file, "a") as fin:
							SeqIO.write(record,fin,"fasta")
							prot_count+=1
				print("".join(["Found ", str(prot_count), " new proteins for ",genome]))

			else:
				print("".join(['miniprot.py did not run correctly for ', genome]))
	else:
		print("".join(['mini-finish.py has already been run please delete or move ', new_orthos]))
else:
	print("".join(['mini-finish.py has already been run please delete or move ', gff_dir]))

##Go through overlapping files
if not os.path.isfile(overlap_out):
	overlap_orthos={}
	for genome in os.listdir(minidir):
		genome_dir=f"{minidir}/{genome}"
		overlap_file=f"{genome_dir}/{genome}_overlapping-orthologs.csv"
		overlap_list=[]
		if os.path.isfile(overlap_file):
			with open(overlap_file,"r") as fin:
				for line in fin:
					orthologs=line.split("\n")[0]
					ortho1=orthologs.split(",")[0]
					ortho2=orthologs.split(",")[1]
					if ortho1 != "OG0000000":
						ortho1_num=int((ortho1.split("OG")[1]).lstrip("0"))
					else:
						ortho1_num=0
					if ortho2 != "OG0000000":
						ortho2_num=int((ortho2.split("OG")[1]).lstrip("0"))
					else:
						ortho2_num=0

					if ortho1_num > ortho2_num:
						ortho_rearranged=f"{ortho2},{ortho1}"
					if ortho1_num < ortho2_num:
						ortho_rearranged=f"{ortho1},{ortho2}"
					if ortho_rearranged not in overlap_list:
						overlap_list.append(ortho_rearranged)
						if ortho_rearranged not in overlap_orthos:
							overlap_orthos[ortho_rearranged]=0
						overlap_orthos[ortho_rearranged]+=1


	#Write ortholog file count
	if len(overlap_orthos)>1:
		with open(overlap_out,"w") as f:
			f.write("".join(['ortholog_1,ortholog2,count\n']))
			for key in overlap_orthos:
				f.write("".join([key,",",str(overlap_orthos[key]),"\n"]))
