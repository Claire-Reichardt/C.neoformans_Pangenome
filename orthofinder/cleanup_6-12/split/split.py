##Last updated:	5/13/26

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
parser = argparse.ArgumentParser(description='dirty script for separating combined proteins by CNAG proteins on GACRC Slurm',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('working_directory', help="directory containing: directories of orthologs annotated as accessories with both orthogroup file inside and a CNAGs.fa file of CNAG proteins to split orthogroups on")
parser.add_argument('--submit', help='Submit BLAST job directly to slurm', action='store_true')
args = parser.parse_args()

##file definitions
userdir=os.getcwd()
wrk_dir=os.path.abspath(args.working_directory)

##Write shell to blast combo ortholog
for orthogroup in os.listdir(wrk_dir):
	if os.path.isdir(orthogroup):
		orthodir=f"{wrk_dir}/{orthogroup}/"
		CNAGs=f"{orthodir}CNAGs.fa"
		split_cnag=None
		split_ortho=None
		acc_cnag=None
		acc_ortho=None
		for record in SeqIO.parse(CNAGs,"fasta"):
			ortho_id=(record.description).split("\t")[1]
			if ortho_id!=orthogroup:
				split_cnag=record.id
				split_ortho=ortho_id
			else:
				acc_cnag=record.id
				acc_ortho=ortho_id

		#Now make blast shell file
		blast_dir=f"{orthodir}blastp/"
		sub_file=f"{orthodir}{split_ortho}.fa"
		os.makedirs(os.path.dirname(blast_dir), exist_ok=True)
		blastout=f"{orthodir}{orthogroup}_blastp.out"
		if not os.path.isfile(blastout):
			shellfile=f"{blast_dir}blastp.sh"
			os.makedirs(os.path.dirname(shellfile),exist_ok=True)
			shell=open(shellfile,"w")
			shell.write("""#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=""")
			shell.write(blast_dir)
			shell.write("""%x_%j.out
""")
			shell.write('\n')
			shell.write('\n')

			#write load modules
			shell.write("""module purge
ml Python/3.10.4-GCCcore-11.3.0
ml BLAST+/2.13.0-gompi-2022a""")
			shell.write("".join(['\n', '\n']))

			#change directories
			shell.write("".join(['cd ', orthodir, '\n', '\n']))

			#run blast
			shell.write("".join(['makeblastdb -in ', sub_file,' -dbtype prot -out ', blast_dir,
			split_ortho, '\n']))

			shell.write("".join(['blastp -query ', CNAGs, ' -db ', blast_dir, split_ortho,
			' -out ',blastout,
			' -outfmt "6 qseqid sseqid sseq evalue" -evalue 1e-05', '\n']))

			shell.close()
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile, '\n']))
				os.system(execute)

		##go through output
		if os.path.isfile(blastout):
			split_og=f"{orthodir}{split_ortho}.fa"
			split_og_ids=[]
			split_og_seqs=[]
			missing_og=f"{orthodir}{acc_ortho}.fa"
			missing_og_ids=[]
			missing_og_seqs=[]
			if os.path.isfile(missing_og):
				for record in SeqIO.parse(missing_og,"fasta"):
					missing_og_ids.append(record.id)
					missing_og_seqs.append(record)

			with open(blastout,'r') as fin:
				for line in fin:
					hit=line.split("\t")[0]
					id=line.split("\t")[1]
					seq=Seq(line.split("\t")[2])
					if "-" in seq:
						seq=Seq(seq.replace("-",""))
					new_id=f"{id}_B"
					if hit==acc_cnag:
						if id not in missing_og_ids:
							if new_id not in missing_og_ids:
								new_record=SeqRecord(seq,id=new_id,name="",description="")
								missing_og_seqs.append(new_record)
								missing_og_ids.append(new_id)
			with open(blastout, 'r') as fin:
				for line in fin:
					hit=line.split("\t")[0]
					id=line.split("\t")[1]
					seq=Seq(line.split("\t")[2])
					if "-" in seq:
						seq=Seq(seq.replace("-",""))
					new_id=f"{id}_B"
					second_id=f"{id}_A"
					if hit==split_cnag:
						if new_id in missing_og_ids:
							if second_id not in split_og_ids:
								new_record=SeqRecord(seq,id=second_id,name="",description="")
								split_og_ids.append(second_id)
								split_og_seqs.append(new_record)


			if os.path.isfile(split_og):
				for record in SeqIO.parse(split_og,"fasta"):
					split_id=f"{record.id}_A"
					if split_id not in split_og_ids:
						split_og_ids.append(record.id)
						split_og_seqs.append(record)

			##now write final record files
			new_split=f"{orthodir}{split_ortho}_new.fa"
			new_acc=f"{orthodir}{acc_ortho}_new.fa"
			SeqIO.write(split_og_seqs,new_split,'fasta')
			SeqIO.write(missing_og_seqs,new_acc,'fasta')
