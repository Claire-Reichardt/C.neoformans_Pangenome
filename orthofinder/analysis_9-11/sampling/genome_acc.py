##Last updated:	10/29/25

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import random
import pandas as pd
parser = argparse.ArgumentParser(description='Code to count isolates with given accessory genes on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
parser.add_argument('classify', metavar='ORTHOLOG_CLASSIFICATION', help='classify.txt from classify.py output')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
corefile=f"{userdir}/core_pergenome_knownclass.tsv"
noncorefile = f"{userdir}/noncore_pergenome_knownclass.tsv"
accfile = f"{userdir}/accessoryonly_pergenome_knownclass.tsv"
orderfile=f"{userdir}/sampling_order_knownclass.tsv"

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

CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.txt is missing or not readable!!")

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

#Get the classification for each ortholog
with open(classfile,'r') as fin:
	lines=fin.readlines()

class_dict={}
for line in lines:
	no_n=line.split("\n")[0]
	ortho=line[0:9]
	gene_class=no_n.split(",")[1]
	class_dict[ortho]=gene_class

##Count number of noncore/accessory orthos when randomly going through files
#create dataframes
core_df = pd.DataFrame()
noncore_df = pd.DataFrame()
acc_df = pd.DataFrame()
order_df=pd.DataFrame()

#iterate 100 times
sample=0

for i in range(1):
	random.shuffle(file_list)
	sample+=1
	core_count=0
	noncore_count=0
	accessory_count=0
	unique_noncore=[]
	unique_accessory=[]
	unique_core=[]
	core_list=[]
	acc_list=[]
	noncore_list=[]
	genome_order=[]
	accessory_count={}

	for file in file_list:
		#get proper genome name
		genomename= file.split("_")[0]
		name = os.path.splitext(genomename)[0]
		genome_order.append(name)
		filename=f"{genomesdir}/{file}"
		acc_list=[]

		#read every line in the file
		with open(filename, 'r') as fin:
			lines=fin.readlines()

		for line in lines:
			ortho=line[0:9]
			if ortho in class_dict:
				if class_dict[ortho]=="singleton":
					if ortho not in unique_noncore:
						unique_noncore.append(ortho)
						noncore_count += 1

				if class_dict[ortho]=="accessory":
					if ortho not in unique_accessory:
						unique_noncore.append(ortho)
						noncore_count+=1
						unique_accessory.append(ortho)
						accessory_count[ortho]=0
					if ortho not in acc_list:
						accessory_count[ortho]+=1
						acc_list.append(ortho)

				if class_dict[ortho]=="core":
					if ortho not in unique_core:
						unique_core.append(ortho)
						core_count+=1

		#add counts for each file
		noncore_list.append(str(noncore_count))
		acc_list.append(str(accessory_count))
		core_list.append(str(core_count))

	#create a column for every loop
	samplecount=str(sample)
	samplenum=f"Sample_{samplecount}"
	core_df[samplenum]= pd.Series(core_list)
	noncore_df[samplenum]= pd.Series(noncore_list)
	acc_df[samplenum]= pd.Series(acc_list)
	order_df[samplenum]=pd.Series(genome_order)

	#close things
	fin.close()

for ortho in accessory_count:
	print("".join([ortho, '\t',str(accessory_count[ortho])]))
#export dataframes to tsv
noncore_df.to_csv(noncorefile, sep='\t')
acc_df.to_csv(accfile, sep='\t')
order_df.to_csv(orderfile, sep='\t')
core_df.to_csv(corefile, sep='\t')
