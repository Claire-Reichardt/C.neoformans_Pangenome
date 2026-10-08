##Last updated:	4/29/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import pandas as pd
parser = argparse.ArgumentParser(description='Code to sort Ortholog groups into core/accessory/singleton on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outfile = f"{userdir}/genefamily_pergenome.tsv"
classfile= f"{userdir}/classify.csv"
countfile=f"{userdir}/ortho_count.csv"

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
countentries=0
for file in os.listdir(genomesdir):
	countentries +=1
total_genomes=str(countentries)
print("".join(['Found ', total_genomes, ' genomes']))

#define cutoffs
core = round(int(total_genomes) * 0.95)
print("".join(['The core gene count cutoff is: ', str(core)]))

singleton = 1

#make a dictionary of orthologues and their gene classifcation based on cutoffs
ortho_counter={}

with open(classfile, 'a') as f:

	for file in os.listdir(genomesdir):
		unique_orthos=[]
		genomename= file.split("_")[0]
		name = os.path.splitext(genomename)[0]
		filename=f"{genomesdir}/{file}"

		#read every line in the file
		with open(filename, 'r') as fin:
			lines=fin.readlines()

		for line in lines:
			ortho=line[0:9]
			if ortho not in unique_orthos:
				unique_orthos.append(ortho)
				if ortho not in ortho_counter:
					ortho_counter[ortho]=0
				ortho_counter[ortho] +=1

#classify based on counts
dict_orthos={}

for ortho in ortho_counter:
	if ortho_counter[ortho]<= core:
		dict_orthos[ortho]="accessory"
	elif ortho_counter[ortho]>core:
		dict_orthos[ortho]="core"

#print to files
dict_orthos_df=pd.DataFrame.from_dict(data=dict_orthos, orient='index')
dict_orthos_df.to_csv(classfile, sep=',')

ortho_counter_df=pd.DataFrame.from_dict(data=ortho_counter, orient='index')
ortho_counter_df.to_csv(countfile, sep=',')

#Count number of core/accessory/singleton genes in each genome
with open(outfile, 'a') as f:
	f.write("".join(['name','\t','core_count','\t','accessory_count','\n']))

for file in os.listdir(genomesdir):
	core_count=0
	acc_count=0

	genomename= file.split("_")[0]
	name = os.path.splitext(genomename)[0]
	filename=f"{genomesdir}/{file}"

	if "CNAG" not in name:
		with open(filename, 'r') as fin:
			lines=fin.readlines()
		for line in lines:
			ortho=line[0:9]
			if ortho in dict_orthos:
				geneclass = dict_orthos[ortho]
				if geneclass == "accessory":
					acc_count += 1

				elif geneclass == "core":
					core_count += 1
	if "CNAG" in name:
		with open(filename, 'r') as fin:
			lines=fin.readlines()
		for line in lines:
			ortho=line[0:9]
			if ortho in dict_orthos:
				geneclass = dict_orthos[ortho]
				if geneclass == "accessory":
					acc_count += 1

				elif geneclass == "core":
					core_count += 1

	with open(outfile, 'a') as f:
		f.write("".join([name,'\t',str(core_count),'\t', str(acc_count),'\n']))

fin.close()
f.close()
