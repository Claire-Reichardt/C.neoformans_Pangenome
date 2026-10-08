##Last updated:	7/23/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to randomly sample genomes and compare their orthologs to singular other genomes on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
compfile = f"{userdir}/comparison.tsv"

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
count = 0
for fp in os.listdir(orthodir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any Orthogroup files in: ', orthodir, '. Please check OrthoFinder ran completely in provided directory.']))
else:
	print("".join(['Found ', str(count), ' Orthogroup files in: ', orthodir]))

#Get total genomes used for analysis
file_list= []
countentries=0
for file in os.listdir(genomesdir):
	countentries +=1
	file_list.append(file)
total_genomes=str(countentries)
print("".join(['Found ', str(countentries), ' genomes in: ', genomesdir]))

##Count differences in orthos when randomly going through files
#create dataframes
compare_df=pd.DataFrame()

for file in file_list:
	#make lists
	unique_orthos_A=[]
	percent_shared=[]

	#get proper genome name
	genomename_A= file.split("_")[0]
	name_A = os.path.splitext(genomename_A)[0]
	filename_A=f"{genomesdir}/{file}"

	#read every line in the file
	with open(filename_A, 'r') as fin:
		lines=fin.readlines()

	for line in lines:
		ortho=line[0:9]
		if ortho not in unique_orthos_A:
			unique_orthos_A.append(ortho)

	#open every other file
	for file in file_list:
		unique_orthos_B=[]
		#similarity_count=0
		genomename_B= file.split("_")[0]
		name_B = os.path.splitext(genomename_B)[0]
		filename_B=f"{genomesdir}/{file}"

		with open(filename_B, 'r') as fin_B:
			lines_B=fin_B.readlines()
		for line in lines_B:
			ortho_B=line[0:9]
			if ortho_B not in unique_orthos_B:
				unique_orthos_B.append(ortho_B)

		shared=set(unique_orthos_A)& set(unique_orthos_B)
		similarity_count=len(shared)
		percent=(similarity_count*2 / (len(unique_orthos_A)+len(unique_orthos_B)))*100
		percent_shared.append(str(percent))

	#create a column for every loop
	compare_df[name_A]= pd.Series(percent_shared)

#export dataframes to tsv
compare_df.to_csv(compfile, sep='\t')
