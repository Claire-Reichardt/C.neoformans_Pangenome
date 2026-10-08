##Last updated:	2/6/25

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
parser = argparse.ArgumentParser(description='Script for separating combined proteins by CNAG proteins on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit BLAST job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
required.add_argument('-p',dest='proteindir', metavar='PROTEIN_DIRECTORY', help='Protein files (.protein.fa) folder from Funannotate output given to OrthoFinder as input', required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
separatedir=f"{userdir}/separate_byCNAGs/"

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

PROTEINDIR = args.proteindir
if os.path.isdir(PROTEINDIR) and os.access(PROTEINDIR, os.R_OK):
	proteindir = (os.path.abspath(args.proteindir))
else:
	sys.exit("ERROR: Protein sequence directory is missing or not readable!!")

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


##Go through and find orthologs that have two representatives for some but one for others
blast_orthos=[]
blast_orthos_CNAGs={}

if not os.path.isdir(separatedir):
	for file in os.listdir(orthodir):
		filename=("".join([orthodir, '/', file]))
		orthogroup = os.path.splitext(file)[0]
		CNAGs=[]
		for record in SeqIO.parse(filename, "fasta"):
			if "CNAG" in record.id:
				CNAGs.append(record.id)
		if len(CNAGs)>1:
			for CNAG1 in CNAGs:
				num1=int(((CNAG1.split("-t")[0]).split("_")[1]).lstrip("0"))
				for CNAG2 in CNAGs:
					if CNAG1 != CNAG2:
						num2=int(((CNAG2.split("-t")[0]).split("_")[1]).lstrip("0"))
						if num2 ==(num1+1) or num2==(num1-1):
							if orthogroup not in blast_orthos:
								blast_orthos.append(orthogroup)
								CNAG_list=[CNAG1,CNAG2]
								blast_orthos_CNAGs[orthogroup]=CNAG_list

	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found ',str(len(blast_orthos_CNAGs.keys())), ' orthologs to split.']))
else:
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found existing directory: ', separatedir, ' Continuing.']))

#Now set up blasts for all sequences
if not os.path.isdir(separatedir):
	os.makedirs(os.path.dirname(separatedir),exist_ok=True)
	for orthogroup in blast_orthos:
		protein_list=list(blast_orthos_CNAGs[orthogroup])
		blast_seqs=f"{separatedir}{orthogroup}/{orthogroup}_queries.fa"
		protein_seqs=[]
		if not os.path.isfile(blast_seqs):
			filename=f"{orthodir}/{orthogroup}.fa"
			if os.path.isfile(filename):
				for record in SeqIO.parse(filename, "fasta"):
					if record.id in protein_list:
						protein_seqs.append(record)
			os.makedirs(os.path.dirname(blast_seqs),exist_ok=True)
			SeqIO.write(protein_seqs,blast_seqs,'fasta')

	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Created ',str(len(blast_orthos_CNAGs.keys())),
	' CNAG protein query files to check for similarities.']))
else:
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found existing directory: ', separatedir, ' Continuing.']))

for orthogroup in os.listdir(separatedir):
	subjectfile=f"{separatedir}{orthogroup}/{orthogroup}_queries.fa"
	seqfile=f"{orthodir}/{orthogroup}.fa"
	database=f"{separatedir}{orthogroup}/db/"
	os.makedirs(os.path.dirname(database), exist_ok=True)
	blastout=f"{separatedir}{orthogroup}/{orthogroup}_blastp.out"
	if not os.path.isfile(blastout):
		if os.path.isfile(subjectfile) and os.path.isfile(seqfile):
			shellfile=f"{separatedir}{orthogroup}/blastshell.sh"
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
			shell.write(f"{separatedir}{orthogroup}/")
			shell.write("""%x_%j.out
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
			shell.write("".join(['cd ', separatedir, orthogroup, '\n', '\n']))

			#run blast
			shell.write("".join(['makeblastdb -in ', subjectfile,' -dbtype prot -out ', database,
			orthogroup, '\n']))

			shell.write("".join(['blastp -query ', seqfile, ' -db ', database, orthogroup,
			' -out ',blastout,
			' -outfmt "6 qseqid sseqid qstart qend qcovs gapopen evalue bitscore" -evalue 1e-05', '\n']))

			shell.close()
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile, '\n']))
				os.system(execute)

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Completed BLASTs. Continuing.']))

#Now check if proteins align to both representatives for splitting
split_proteins={}
orthogroups_split=[]
for orthogroup in os.listdir(separatedir):
	blastout=f"{separatedir}{orthogroup}/{orthogroup}_blastp.out"
	if os.path.isfile(blastout):
		protein_count={}
		with open(blastout,'r') as fin:
			lines=fin.readlines()
		for line in lines:
			qseqid=line.split("\t")[0]
			sseqid=line.split("\t")[1]
			qstart=line.split("\t")[2]
			qend=line.split("\t")[3]
			qcovs=line.split("\t")[4]
			gapopen=line.split("\t")[5]
			evalue=line.split("\t")[6]
			bitscore=(line.split("\t")[7]).split("\n")[0]
			if int(qcovs) < 90:
				if Decimal(evalue) < Decimal('1.0e-50'):
					if qseqid not in protein_count:
						protein_count[qseqid]=0
					protein_count[qseqid]+=1
		for protein in protein_count:
			split_sites=[]
			split_lines=[]
			if protein_count[protein]>1:
				groups=[]
				for line in lines:
					if protein in line:
						sseqid=line.split("\t")[1]
						qstart=line.split("\t")[2]
						qend=line.split("\t")[3]
						gapopen=line.split("\t")[5]
						groups.append(sseqid)

				for group in groups:
					for line in lines:
						if group and protein in line:
							qstart=line.split("\t")[2]
							qend=line.split("\t")[3]
							evalue=line.split("\t")[6]
							if int(qend)-int(qstart)>=35:
								if Decimal(evalue) < Decimal('1.0e-50'):
									split_site=f"{qstart}:{qend}"
									if split_site not in split_sites:
										split_sites.append(split_site)
										split_lines.append(line)
			if len(split_sites)>1:
				appended_split_sites=[]
				for split_1 in split_sites:
					start_1=int(split_1.split(":")[0])
					end_1=int(split_1.split(":")[1])
					for split_2 in split_sites:
						start_2=int(split_2.split(":")[0])
						end_2=int(split_2.split(":")[1])
						if split_1 !=split_2:
							if end_1 > start_2 and start_1 < start_2:
								if end_2 > end_1+35:
									for line_1 in split_lines:
										gapopen_1=int(line_1.split("\t")[5])
									for line_2 in split_lines:
										gapopen_2=int(line_2.split("\t")[5])
										if line_1 !=line_2:
											if gapopen_1 > gapopen_2:
												start_2a=(start_2)-1
												append_split_1=f"{start_1}:{start_2a}"
												append_split_2=f"{start_2}:{end_2}"
												if append_split_1 not in appended_split_sites:
													appended_split_sites.append(append_split_1)
												if append_split_2 not in appended_split_sites:
													appended_split_sites.append(append_split_2)
											if gapopen_1 < gapopen_2:
												end_1b=(end_1)-1
												append_split_1=f"{start_1}:{end_1}"
												append_split_2=f"{end_1b}:{end_2}"
												if append_split_1 not in appended_split_sites:
													appended_split_sites.append(append_split_1)
												if append_split_2 not in appended_split_sites:
													appended_split_sites.append(append_split_2)


							if end_1 < start_2:
								if split_1 not in appended_split_sites:
									appended_split_sites.append(split_1)
								if split_2 not in appended_split_sites:
									appended_split_sites.append(split_2)
							if end_1==start_2:
								if split_1 not in appended_split_sites:
									appended_split_sites.append(split_1)
								end_1b=(end_1)-1
								append_split_2=f"{end_1b}:{end_2}"
								if append_split_2 not in appended_split_sites:
									appended_split_sites.append(append_split_2)

				if len(appended_split_sites)>1:
					split_proteins[protein]=appended_split_sites
					if orthogroup not in orthogroups_split:
						orthogroups_split.append(orthogroup)
now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Found ',str(len(orthogroups_split)),
' orthogroups to to split based on protein alignments. Continuing. ']))

#write new split fastas into a new file
all_split_count=0
numbers = [1,2,3,4,5]
alphabet = string.ascii_uppercase
mapping = {number: alphabet[i] for i, number in enumerate(numbers)}
new_prots_dir=f"{userdir}/split_protein_files/"
if not os.path.isdir(new_prots_dir):
	os.makedirs(os.path.dirname(new_prots_dir),exist_ok=True)
	for file in os.listdir(proteindir):
		filename=("".join([proteindir, '/', file]))
		if "cryptococcus_" in file:
			genome=(file.split("cryptococcus_")[1]).split(".proteins.fa")[0]
		else:
			genome=file.split(".fa")[0]
		out_filename=f"{new_prots_dir}{genome}.fa"
		protein_seqs=[]

		for record in SeqIO.parse(filename,"fasta"):
			if record.id in split_proteins:
				cuts=list(split_proteins[record.id])
				cut_count=0
				for cut in cuts:
					cut_count+=1
					all_split_count+=1
					letter=mapping[cut_count]
					id=f"{record.id}_{letter}"
					start=int(cut.split(":")[0])
					end=int(cut.split(":")[1])
					cut_seq=record.seq[start:end]
					new_record=SeqRecord(cut_seq,id,"",record.id)
					protein_seqs.append(new_record)
			else:
				protein_seqs.append(record)

		SeqIO.write(protein_seqs, out_filename,'fasta')

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Split ', str(all_split_count), ' proteins. New protein files located in: ', new_prots_dir]))

#Now check if all orthologs have been split accordingly

for file in os.listdir(orthodir):
	filename=("".join([orthodir, '/', file]))
	orthogroup = os.path.splitext(file)[0]
	if orthogroup in orthogroups_split:
		blast_seqs=f"{separatedir}{orthogroup}/{orthogroup}_queries.fa"
		CNAG_combo_length=0
		for record in SeqIO.parse(blast_seqs,"fasta"):
			if "CNAG" in record.id:
				CNAG_combo_length+=(len(record.seq))
		for record_og in SeqIO.parse(filename,"fasta"):
			if len(record_og.seq) >= CNAG_combo_length*0.75:
				header=(record_og.id).split("_")[0]
				if len(header)<10:
					genome=header
					new_filename=f"{new_prots_dir}{genome}.fa"
				if len(header)>=10:
					for file_2 in os.listdir(new_prots_dir):
						genome=file_2.split(".fa")[0]
						if header in genome:
							new_filename=f"{new_prots_dir}{genome}.fa"
				protein_count=0
				for record_new in SeqIO.parse(new_filename,"fasta"):
					if record_og.id in record_new.id:
						protein_count+=1
				if protein_count<2:
					print(record_og.id)
