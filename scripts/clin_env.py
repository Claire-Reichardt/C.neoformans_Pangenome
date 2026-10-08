##Last updated:	8/14/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description='Code to sort Ortholog groups presence in clinical or environmental strains on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
parser.add_argument('envkey', metavar='ENVIRONMENT_KEY', help='.CSV file containing header used for proteins of each genome and the clinical/environmental isolation associated with it')
parser.add_argument('classify', metavar='ORTHOLOG_CLASSIFICATION', help='classify.txt from classify.py output')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outfile = f"{userdir}/clinvsenv.csv"

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
GENOMESDIR =f"{args.orthodir}/all_orthos"
ENVKEY=args.envkey
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")
if os.path.isfile(ENVKEY) and os.access(ENVKEY, os.R_OK):
	envkey = (os.path.abspath(ENVKEY))
else:
	sys.exit("ERROR: matkey.csv provided is missing or not readable!!")

CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.txt is missing or not readable!!")

#Get the classification for each ortholog
with open(classfile,'r') as fin:
	lines=fin.readlines()

class_dict={}
for line in lines:
	no_n=line.split("\n")[0]
	ortho=line[0:9]
	gene_class=no_n.split(",")[1]
	class_dict[ortho]=gene_class

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
	if "CNAG" not in file:
		countentries +=1
total_genomes=str(countentries)

#Read input csv to get key
env_dict={}
with open(envkey, 'r') as f:
	for line in f:
		items=line.split(',')
		key, values = items[0], items[1]
		values_2= values.split('\n')[0]
		env_dict[key]=values_2

#Count number of core/accessory/singleton genes in each genome
with open(outfile, 'a') as f:
	f.write("".join(['orthogroup',',','classification',',','clinical_count', ',','environmental_count', '\n']))

for file in os.listdir(orthodir):
	clinical_count=0
	environmental_count=0
	unique_genome=[]

	filename=("".join([orthodir, '/', file]))
	orthogroup = os.path.splitext(file)[0]

	if orthogroup in class_dict:
		geneclass = class_dict[orthogroup]

		for record in SeqIO.parse(filename, "fasta"):
			name=record.id
			base=name.split('_')[0]
			if base in env_dict:
				isolation= env_dict[base]
				if base not in unique_genome:
					if isolation =="Clinical":
						unique_genome.append(base)
						clinical_count +=1
					if isolation=="Environmental":
						unique_genome.append(base)
						environmental_count +=1

	with open(outfile, 'a') as f:
		f.write("".join([orthogroup,',', geneclass ,',', str(clinical_count), ',', str(environmental_count),'\n']))

