##Last updated:	5/14/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to count number of proteins per genome in each orthogroup',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
GENOMESDIR =f"{args.orthodir}/all_orthos"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")

##print files
ortho_count = []
for fp in os.listdir(orthodir):
	ext = os.path.splitext(fp)[-1].lower()
	ortho = os.path.splitext(fp)[0]
	if ext == ".fa":
		ortho_count.append(ortho)
if ortho_count ==0:
	sys.exit("".join(['ERROR: Unable to find any Orthogroup files in: ', orthodir, '. Please check OrthoFinder ran completely in provided directory.']))
else:
	print("".join(['Found ', str(len(ortho_count)), ' Orthogroup files in: ', orthodir]))

ortho_count.sort()

#Get total genomes used for analysis
file_list= []
countentries=0
for file in os.listdir(genomesdir):
	countentries +=1
	file_list.append(file)
total_genomes=str(countentries)
print("".join(['Found ', str(countentries), ' genomes in: ', genomesdir]))

##file definitions
userdir = os.getcwd()
compfile = f"{userdir}/copy_num.csv"

#create dataframes
copy_df=pd.DataFrame()
copy_df["OGs"]=pd.Series(ortho_count)

for file in file_list:
	final_count=[]
	prot_count={}
	name= os.path.splitext(file)[0]
	filename=f"{genomesdir}/{file}"
	with open(filename,'r') as fin:
		lines=fin.readlines()
	for line in lines:
		ortho=line[0:9]
		if ortho not in prot_count:
			prot_count[ortho]=0
		prot_count[ortho]+=1
	for ortho in ortho_count:
		if ortho in prot_count:
			final_count.append(prot_count[ortho])
		else:
			final_count.append(0)
	copy_df[name]=pd.Series(final_count)

copy_df.to_csv(compfile,sep=',')
