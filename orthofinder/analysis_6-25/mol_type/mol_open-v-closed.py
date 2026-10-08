##Last updated:	8/17/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to randomly sample Ortholog groups per molecular type-specific genomes and classify them per genome GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
parser.add_argument('molkey', metavar='MOLECULARTYPE_KEY', help='.TSV file containing strain name and the molecular type associated with it')
parser.add_argument('moltype',metavar='MOLECULARTYPE_SEARCH', nargs=1,help='name of molecular type to randomly sample')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
newortho = f"{userdir}/{str(args.moltype[0])}_newgene.csv"

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
GENOMESDIR =f"{args.orthodir}/all_orthos"
MOLKEY=args.molkey
if os.path.isfile(MOLKEY) and os.access(MOLKEY, os.R_OK):
	molkey = (os.path.abspath(MOLKEY))
else:
	sys.exit("ERROR: molkey.tsv provided is missing or not readable!!")
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")

#Read input csv to get key
clade_dict={}
with open(molkey, 'r') as f:
	for line in f:
		items=line.split('\t')
		key, values = items[0], items[1]
		values_2= values.split('\n')[0]
		clade_dict[key]=values_2

#Get total genomes used for analysis
file_list= []
countentries=0
for file in os.listdir(genomesdir):
	genomename=file.split("_")[0]
	name=os.path.splitext(genomename)[0]
	if "CNAG" not in file:
		if name in clade_dict:
			if clade_dict[name]==str(args.moltype[0]):
				countentries +=1
				file_list.append(file)
total_genomes=str(countentries)
print("".join(['Found ', str(countentries), ' genomes that are ', str(args.moltype[0]),' in: ', genomesdir]))

##Count number of noncore/accessory orthos when randomly going through files
#create dataframe
gene_df=pd.DataFrame()

#iterate 100 times
sample=0

for i in range(100):
	random.shuffle(file_list)
	sample+=1
	ortho_counter={}
	new_genes=[]

	for file in file_list:
		#get proper genome name
		if "_" in file:
			genomename= file.split("_")[0]
		else:
			genomename=file
		name = os.path.splitext(genomename)[0]
		filename=f"{genomesdir}/{file}"

		#read every line in the file
		with open(filename, 'r') as fin:
			lines=fin.readlines()

		for line in lines:
			ortho=line[0:9]
			if ortho in ortho_counter:
				ortho_counter[ortho].append(name)
			else:
				ortho_counter[ortho]=[name]

		new_genes.append(len(ortho_counter))

	#create a column for every loop
	samplecount=str(sample)
	samplenum=f"Sample_{samplecount}"
	gene_df[samplenum]= pd.Series(new_genes)

	#close things
	fin.close()

#export dataframes to csv
gene_df.to_csv(newortho, sep=',')
