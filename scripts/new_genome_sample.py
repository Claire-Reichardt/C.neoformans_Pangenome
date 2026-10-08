##Last updated:	6/5/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to randomly sample Ortholog groups per genome and classify them per genome GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
noncorefile = f"{userdir}/noncore_pergenome.csv"
accfile = f"{userdir}/accessory_pergenome.csv"

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

#Get total genomes used for analysis
file_list= []
countentries=0
for file in os.listdir(genomesdir):
	countentries +=1
	file_list.append(file)
total_genomes=str(countentries)
print("".join(['Found ', str(countentries), ' genomes in: ', genomesdir]))

##Count number of noncore/accessory orthos when randomly going through files
#create dataframes
noncore_df = pd.DataFrame()
acc_df = pd.DataFrame()

#iterate 100 times
sample=0

for i in range(100):
	random.shuffle(file_list)
	sample+=1
	sing_count=0
	accessory_count=0
	unique_sing=[]
	unique_accessory=[]
	unique_core=[]
	noncore_list=[]
	acc_list=[]
	genomes=0
	ortho_counter={}

	for file in file_list:
		#define cutoffs
		genomes+=1
		core=round(genomes * 0.95)
		singleton=1

		#get proper genome name
		genomename= file.split("_")[0]
		name = os.path.splitext(genomename)[0]
		filename=f"{genomesdir}/{file}"

		#read every line in the file
		with open(filename, 'r') as fin:
			lines=fin.readlines()

		for line in lines:
			ortho=line[0:9]
			if ortho not in ortho_counter:
				ortho_counter[ortho]=0
			ortho_counter[ortho] +=1

			if ortho in ortho_counter:
				if singleton < ortho_counter[ortho] < core:
					if ortho not in unique_accessory:
						unique_accessory.append(ortho)
						accessory_count+=1
					if ortho in unique_sing:
						unique_sing.remove(ortho)
						if sing_count>0:
							sing_count-=1
					if ortho in unique_core:
						unique_core.remove(ortho)

				if ortho_counter[ortho]==singleton:
					if ortho not in unique_sing:
						unique_sing.append(ortho)
						sing_count+=1

				if ortho_counter[ortho]>=core:
					if ortho not in unique_core:
						unique_core.append(ortho)
					if ortho in unique_sing:
						unique_sing.remove(ortho)
						if sing_count>0:
							sing_count-=1
					if ortho in unique_accessory:
						unique_accessory.remove(ortho)
						if accessory_count>0:
							accessory_count-=1

		#add counts for each file
		noncore_list.append(str(accessory_count+sing_count))
		acc_list.append(str(accessory_count))

	#create a column for every loop
	samplecount=str(sample)
	samplenum=f"Sample_{samplecount}"
	noncore_df[samplenum]= pd.Series(noncore_list)
	acc_df[samplenum]= pd.Series(acc_list)

	#close things
	fin.close()

#export dataframes to csv
noncore_df.to_csv(noncorefile, sep=',')
acc_df.to_csv(accfile, sep=',')
