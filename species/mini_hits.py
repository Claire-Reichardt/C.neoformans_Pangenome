##Last updated: 8/25/26

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
import pandas as pd
parser = argparse.ArgumentParser(description='Submodule of species_mini.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('mini_dir')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
mini_dir=os.path.abspath(args.mini_dir)
outfile=f"{userdir}/species_miniprot_identities.csv"

##Go through all directories
genome_hit_dict={}
for genome in os.listdir(mini_dir):
	best_orthos={}
	genome_dir=f"{mini_dir}/{genome}"
	for file in os.listdir(genome_dir):
		if ".gff3" in file:
			best_hit=0
			file_path=f"{genome_dir}/{file}"
			ortholog=((file.split(genome)[1]).split(".gff3")[0]).split("_")[1]
			with open(file_path,'r') as fin:
				for line in fin:
					if "mRNA" in line:
						att=line.split("\t")[8]
						id=att.split(";")[0]
						id_num=id.split("ID=")[1]
						ortho=id.split("M")[0]
						perid=att.split(";")[2]
						perid_only=float(perid.split("Identity=")[1])
						if perid_only > best_hit:
							best_hit=perid_only
			if best_hit !=0:
				best_orthos[ortholog]=best_hit
	genome_hit_dict[genome]=pd.Series(best_orthos)
df=pd.DataFrame(genome_hit_dict)
df.to_csv(outfile,sep=",")
