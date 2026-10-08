##Last updated:	6/14/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import pandas as pd
parser = argparse.ArgumentParser(description='Code to extract annotation of genes from longest_sequence.py output and from funannotate .gff on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument('gfffile', metavar='GFF3_FILE', help='.gff3 file containing all .gff3 annotations for all genomes')
parser.add_argument('longestseq', metavar='LONGEST_SEQUENCE_FILE', help='.fa file from longest_sequence.py containing representatives for each Ortholog')
args = parser.parse_args()

##file definitions
userdir=os.getcwd()
findfile=f"{userdir}/find_proteins.txt"
outfile=f"{userdir}/ortho_annotations.tsv"
shellfile=f"{userdir}/find_shell.sh"
tmpfindfile=f"{userdir}/find_proteins_partial.txt"

##Check files exist
GFF3_FILE = args.gfffile
if os.path.isfile(GFF3_FILE) and os.access(GFF3_FILE, os.R_OK):
	gff = (os.path.abspath(args.gfffile))
else:
	sys.exit("ERROR: GFF3 file is missing or not readable!")

LONGEST_FILE = args.longestseq
if os.path.isfile(LONGEST_FILE) and os.access(LONGEST_FILE, os.R_OK):
	seq = (os.path.abspath(args.longestseq))
else:
	sys.exit("ERROR: Longest sequence file is missing or not readable!")

#Get ids to search for from fasta file
	#Make sure the GFF3 file has NO spaces-just tabs.
		#I used: sed 's/ /_/g' all_gffs.gff > no_space.gff3

proteinortho_dict={}
records=[]
for record in SeqIO.parse(seq, "fasta"):
	base=record.description
	proteins=base.split(" ")[1]
	records.append(proteins)
	ortholog=record.id
	proteinortho_dict[proteins]=ortholog
sorted_records=sorted(records)
print("".join(['Found: ', str(len(sorted_records)), ' ortholog representative sequences in: ', seq]))

#Make a shell for large searches
if args.submit is True:
	shell=open(shellfile, "w")

	shell.write("""#!/bin/bash
#SBATCH --job-name=find
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
""")

	for protein in sorted_records:
		find=("".join([
				"awk '/",
				protein, "/ &&/mRNA/ {print$9}' ", gff, " >> ", findfile, "\n"]))
		shell.write(find)

#submit shell
	shell.close()
	os.system("".join(['sbatch ', shellfile, '\n']))

#Attach orthologs to annotations
if args.submit is False:
	annots_dict={}

	with open(tmpfindfile, 'r') as fin:
		lines=fin.readlines()
	for line in lines:
		genome=line.split("=")[1]
		protein=genome.split(";")[0]
		if protein in proteinortho_dict:
			ortho=proteinortho_dict[protein]
			annotation=line.split(";",2)[2]
			annot=annotation.split("\n")[0]
			annots_dict[ortho]=annot

	#Export dataframe to tsv
	annots_dict_sorted=dict(sorted(annots_dict.items()))
	annotation_df=pd.DataFrame.from_dict(data=annots_dict_sorted, orient='index')
	annotation_df.to_csv(outfile, sep='\t')
