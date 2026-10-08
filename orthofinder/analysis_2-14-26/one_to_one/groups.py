##Last updated:	3/6/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to compare molecular type orthologs on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)', required=True)
required.add_argument('-k',dest='molkey', metavar='MOLECULARTYPE_KEY', help='.TSV file containing strain name and the molecular type associated with it',required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
compfile = f"{userdir}/mol_comparison.tsv"

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

MOLKEY=args.molkey
if os.path.isfile(MOLKEY) and os.access(MOLKEY, os.R_OK):
	molkey = (os.path.abspath(MOLKEY))
else:
	sys.exit("ERROR: molkey.tsv provided is missing or not readable!!")

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

#Read input csv to get key
clade_dict={}
with open(molkey, 'r') as f:
	for line in f:
		items=line.split('\t')
		key, values = items[0], items[1]
		values_2= values.split('\n')[0]
		clade_dict[key]=values_2

##Count differences in orthos when randomly going through files
#create dataframes
compare_df=pd.DataFrame()

mol_type_df={}

for file in file_list:
	#get proper genome name
	genomename= file.split("_")[0]
	name = os.path.splitext(genomename)[0]
	filename=f"{genomesdir}/{file}"

	if name in clade_dict:
		vn_type=clade_dict[name]

		with open(filename, 'r') as fin:
			lines=fin.readlines()

		for line in lines:
			ortho=line[0:9]
			if vn_type in mol_type_df:
				if ortho not in mol_type_df[vn_type]:
					mol_type_df[vn_type].append(ortho)
			else:
				mol_type_df[vn_type]=[ortho]

	else:
		print("".join(['Could not find ', name,' in ', molkey]))

#export dataframes to tsv
#compare_df.to_csv(compfile, sep='\t')

#Now compare orthologs for groups
per_dict={}
for group in mol_type_df:
	orthogroups=mol_type_df[group]
	for group2 in mol_type_df:
		shared=[]
		orthogroups_2=mol_type_df[group2]
		for ortho in orthogroups:
			if ortho in orthogroups_2:
				if ortho not in shared:
					shared.append(ortho)

		total_group=len(set(orthogroups))+len(set(orthogroups_2))
		per_shared=((len(set(shared))*2)/total_group)*100
		print("".join([group,",", group2,",",str(per_shared)]))
