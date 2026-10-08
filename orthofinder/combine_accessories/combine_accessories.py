##Last updated:	8/12/25

##import system
import sys
import subprocess
import os
import argparse
import numpy as np
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import datetime
import string
from decimal import Decimal
parser = argparse.ArgumentParser(description='Script for combining accessory ortholog proteins with larger orthogroups on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit BLAST job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
required.add_argument('-c',dest='count', metavar='ORTHOLOG_COUNT', help='ortho_count.csv from classify.py output', required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()

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

COUNT = args.count
if os.path.isfile(COUNT) and os.access(COUNT, os.R_OK):
	countfile = (os.path.abspath(COUNT))
else:
	sys.exit("ERROR: ortho_count.csv is missing or not readable!!")

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

#Get the genome count for each ortholog
with open(countfile,'r') as fin:
	lines=fin.readlines()

count_dict={}
for line in lines:
	no_n=line.split("\n")[0]
	ortho=line[0:9]
	gene_count=no_n.split(",")[1]
	count_dict[ortho]=gene_count

#define the longest ortholog
largest_longest = f"{userdir}/longest_large_orthos.fa"
small_longest=f"{userdir}/longest_small_orthos.fa"
if not os.path.isfile(largest_longest) and not os.path.isfile(small_longest):
	longest_large_seq=[]
	longest_small_seq=[]
	for file in os.listdir(orthodir):
		filename=("".join([orthodir, '/', file]))
		orthogroup = os.path.splitext(file)[0]
		longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
		if int(count_dict[orthogroup])>20:
			for record in SeqIO.parse(filename, "fasta"):
				id=record.id
				genome=id.split("_")[0]
				if not longest or len(record) > len(longest):
					if genome not in longest.description:
						longest=record
			if longest:
				rename=SeqRecord(longest.seq, orthogroup, "", longest.id)
				longest_large_seq.append(rename)

		if int(count_dict[orthogroup])<=20:
			for record in SeqIO.parse(filename,"fasta"):
				id=record.id
				genome=id.split("_")[0]
				if not longest or len(record) > len(longest):
					if genome not in longest.description:
						longest=record
			if longest:
				rename=SeqRecord(longest.seq, orthogroup,"",longest.id)
				longest_small_seq.append(rename)

	SeqIO.write(longest_large_seq,largest_longest,"fasta")
	SeqIO.write(longest_small_seq,small_longest,"fasta")

#blast orthologs together
blastdir=f"{userdir}/smallest_blast/"
blastout=f"{blastdir}smallest_blast.tsv"
if not os.path.isfile(blastout):
	dbdir=f"{blastdir}db/"
	os.makedirs(os.path.dirname(dbdir),exist_ok=True)
	shellfile=f"{blastdir}/blastshell.sh"
	os.makedirs(os.path.dirname(shellfile), exist_ok=True)
	shell=open(shellfile,"w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=""")
	shell.write("".join([blastdir,'%x_%j.out','\n']))
	if args.email is not None:
		shell.write("""#SBATCH --mail-user=""")
		shell.write(args.email)
		shell.write("""
#SBATCH --mail-type=FAIL""")
	shell.write("".join(['\n','\n']))

	#write load modules
	shell.write("""module purge
ml BLAST+/2.14.1-gompi-2023a""")
	shell.write("".join(['\n','\n']))

	#change directories
	shell.write("".join(['cd ', blastdir, '\n', '\n']))

	#run blast
	shell.write("".join(['makeblastdb -in ',largest_longest, ' -dbtype prot -out ', dbdir, 'longest_db', '\n',
	'\n']))

	shell.write("".join(['blastp -query ', small_longest, ' -db ', dbdir, 'longest_db -out ',
	blastout, ' -outfmt 6 -evalue 1e-15', '\n', '\n']))

	#submit
	shell.close()
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)

#Go through blast output and define where singletons should go
ortholog_pairs={}
if os.path.isfile(blastout):
	small_evalue={}
	with open(blastout,'r') as fin:
		lines=fin.readlines()
	for line in lines:
		#blast names for outfmt 6
		qseqid=line.split("\t")[0]
		sseqid=line.split("\t")[1]
		pident=line.split("\t")[2]
		#length=line.split("\t")[3]
		#mismatch=line.split("\t")[4]
		#gapopen=line.split("\t")[5]
		#qstart=line.split("\t")[6]
		#qend=line.split("\t")[7]
		#sstart=line.split("\t")[8]
		#ssend=line.split("\t")[9]
		evalue=line.split("\t")[10]
		bitscore=(line.split("\t")[11]).split("\n")[0]
		if int(float(pident))>80:
			small=qseqid
			large=sseqid
			if small not in ortholog_pairs:
				ortholog_pairs[small]=large
				small_evalue[small]=Decimal(evalue)
			else:
				if large!=ortholog_pairs[small]:
					if Decimal(evalue) < Decimal(small_evalue[small]):
						ortholog_pairs[small]=large
						small_evalue[small]=evalue

#Reorder the dictionary
concat_dict={}
if os.path.isfile(blastout):
	for small_ortho in ortholog_pairs:
		small_list=[]
		larger_ortho=ortholog_pairs[small_ortho]
		for key,value in ortholog_pairs.items():
			if value == larger_ortho:
				small_list.append(key)
		concat_dict[larger_ortho]=small_list

#Now combine singletons with the larger orthogroup
if os.path.isfile(blastout):
	cat_files=f"{userdir}/tmp/"
	if not os.path.isdir(cat_files):
		os.makedirs(os.path.dirname(cat_files),exist_ok=True)
		for larger_ortho in concat_dict:
			file_list=[]
			larger_ortho_file=f"{orthodir}/{larger_ortho}.fa"
			file_list.append(larger_ortho_file)
			small_list=concat_dict[larger_ortho]
			for small in small_list:
				small_file=f"{orthodir}/{small}.fa"
				file_list.append(small_file)
			output=f"{cat_files}{larger_ortho}.fa"
			if not os.path.isfile(output):
				with open(output, 'w') as outfile:
					for fname in file_list:
						with open(fname) as infile:
							for line in infile:
								outfile.write(line)

			#remove singleton files that were combined
			for small_ortho in  small_list:
				small_file=f"{orthodir}/{small_ortho}.fa"
				os.remove(small_file)

			#replace old ortholog files with combined ortholog files in tmp dir
			os.rename(output, larger_ortho_file)

print("".join(['Found ', str(len(ortholog_pairs)), ' small orthologs to combine with ',
str(len(concat_dict)), ' larger orthologs.']))
