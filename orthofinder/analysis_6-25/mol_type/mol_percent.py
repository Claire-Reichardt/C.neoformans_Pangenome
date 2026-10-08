##Last updated:	8/11/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to identify orthologs that are specific to one molecular type on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)', required=True)
required.add_argument('-k',dest='molkey', metavar='MOLECULARTYPE_KEY', help='.TSV file containing strain name and the molecular type associated with it',required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
per_file=f"{userdir}/percent_moltype.csv"

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
	genomename= file.split("_")[0]
	name = os.path.splitext(genomename)[0]
	if "CNAG" not in file:
		if name in clade_dict:
			countentries +=1
			file_list.append(file)
		else:
			print("".join(['Unable to find: ', name,' in ', molkey]))
total_genomes=str(countentries)
print("".join(['Found ', str(total_genomes), ' genomes in: ', genomesdir]))

##add orthologs for each file and determine if they are unique
VNI_counter={}
VNII_counter={}
VNBI_counter={}
VNBII_counter={}
total_orthos=[]
for file in file_list:
	genomename= file.split("_")[0]
	name = os.path.splitext(genomename)[0]
	filename=f"{genomesdir}/{file}"
	if name in clade_dict:
		vn_type=clade_dict[name]
		with open(filename, 'r') as fin:
			lines=fin.readlines()
		genome_orthos=[]
		for line in lines:
			ortho=line[0:9]
			if ortho not in genome_orthos:
				genome_orthos.append(ortho)
				if ortho not in total_orthos:
					total_orthos.append(ortho)
				if vn_type=="VNI":
					if ortho not in VNI_counter:
						VNI_counter[ortho]=0
					VNI_counter[ortho]+=1
				if vn_type=="VNII":
					if ortho not in VNII_counter:
						VNII_counter[ortho]=0
					VNII_counter[ortho]+=1
				if vn_type=="VNBI":
					if ortho not in VNBI_counter:
						VNBI_counter[ortho]=0
					VNBI_counter[ortho]+=1
				if vn_type=="VNBII":
					if ortho not in VNBII_counter:
						VNBII_counter[ortho]=0
					VNBII_counter[ortho]+=1

#get percentages
percents={}
for ortho in total_orthos:
	VNI_count=0
	VNII_count=0
	VNBI_count=0
	VNBII_count=0
	mol_list=[]
	if ortho in VNI_counter:
		VNI_count+=VNI_counter[ortho]
		mol_list.append("VNI")
	if ortho in VNII_counter:
		VNII_count+=VNII_counter[ortho]
		mol_list.append("VNII")
	if ortho in VNBI_counter:
		VNBI_count+=VNBI_counter[ortho]
		mol_list.append("VNBI")
	if ortho in VNBII_counter:
		VNBII_count+=VNBII_counter[ortho]
		mol_list.append("VNBII")

	total_count=VNI_count+VNII_count+VNBI_count+VNBII_count
	VNI_per=VNI_count/total_count
	VNII_per=VNII_count/total_count
	VNBI_per=VNBI_count/total_count
	VNBII_per=VNBII_count/total_count
	percents[ortho]=[VNI_per,VNII_per,VNBI_per,VNBII_per,str(mol_list)]

#export into files
per_df=pd.DataFrame.from_dict(data=percents, orient='index',
				columns=['VNI_percent','VNII_percent','VNBI_percent','VNBII_percent','mol_list'])
per_df.to_csv(per_file, sep=',')
