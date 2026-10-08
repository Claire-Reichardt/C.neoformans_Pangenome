##Last updated:	10/3/25

##import system
import sys
import subprocess
import os
import argparse
import pandas as pd
import numpy as np
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import datetime
import string
from decimal import Decimal
parser = argparse.ArgumentParser(description='Script for separating combined proteins on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit BLAST job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
iddir=f"{userdir}/perid/"
os.makedirs(os.path.dirname(iddir),exist_ok=True)

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

##print files
count = 0
for fp in os.listdir(orthodir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any Orthogroup files in: ', orthodir, '. Please check OrthoFinder ran completely in provided directory.']))
else:
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found ', str(count), ' Orthogroup files in: ', orthodir]))

##go through orthologs and set up blast
shelldir=f"{iddir}shells/"
os.makedirs(os.path.dirname(shelldir),exist_ok=True)
dbdir=f"{iddir}dbs/"
os.makedirs(os.path.dirname(dbdir),exist_ok=True)
outdir=f"{iddir}blast_out/"
os.makedirs(os.path.dirname(outdir),exist_ok=True)
outputdir=f"{iddir}shell_out/"
os.makedirs(os.path.dirname(outputdir),exist_ok=True)

for file in os.listdir(orthodir):
	filename=f"{orthodir}/{file}"
	orthogroup=os.path.splitext(file)[0]
	output=f"{outdir}{orthogroup}_out.tsv"
	if not os.path.isfile(output):
		database=f"{dbdir}{orthogroup}.db"
		shellfile=f"{shelldir}{orthogroup}.sh"
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
		shell.write(f"outputdir/{orthogroup}_")
		shell.write("""%j.out
""")
		if args.email is not None:
			shell.write("""#SBATCH --mail-user=""")
			shell.write(args.email)
			shell.write("""
#SBATCH --mail-type=FAIL""")
		shell.write('\n')
		shell.write('\n')

		#write load modules
		shell.write("""module purge
ml Python/3.10.4-GCCcore-11.3.0
ml BLAST+/2.13.0-gompi-2022a""")
		shell.write("".join(['\n', '\n']))

		#change directories
		shell.write("".join(['cd ', outputdir, '\n', '\n']))

		#run blast
		shell.write("".join(['makeblastdb -in ', filename,' -dbtype prot -out ', database,'\n']))

		shell.write("".join(['blastp -query ', filename, ' -db ', database,
		' -out ',output,' -outfmt 6', '\n']))

		shell.close()
		if args.submit is True:
			execute=("".join(['sbatch ', shellfile, '\n']))
			os.system(execute)

#get average percent identities
ortho_pident={}
for file in os.listdir(orthodir):
	orthogroup=os.path.splitext(file)[0]
	output=f"{outdir}{orthogroup}_out.tsv"
	pident_list=[]
	if os.path.isfile(output):
		with open(output,'r') as fin:
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
			#evalue=line.split("\t")[10]
			#bitscore=(line.split("\t")[11]).split("\n")[0]

			if qseqid != sseqid:
				pident_list.append(int(float(pident)))

	if len(pident_list)>0:
		average=sum(pident_list)/len(pident_list)
		ortho_pident[orthogroup]=average
		print("".join([orthogroup,'\t',str(average)]))
	if len(pident_list)==0:
		ortho_pident[orthogroup]="singleton"
		print("".join([orthogroup,'\t',"singleton"]))


##print to file
myKeys = list(ortho_pident.keys())
myKeys.sort()
ortho_pident = {i: ortho_pident[i] for i in myKeys}
pident_outfile=f"{userdir}/orthogroup_percentid.tsv"

ortho_pident_df=pd.DataFrame.from_dict(data=ortho_pident, orient='index')
ortho_pident_df.to_csv(pident_outfile, sep=',')
