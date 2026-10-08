##Last updated:	2/3/25

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
parser = argparse.ArgumentParser(description='Script for separating combined proteins on GACRC Slurm',
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
separatedir=f"{userdir}/separate/"

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
blast_orthos_counts={}

if not os.path.isdir(separatedir):
	for file in os.listdir(orthodir):
		protein_count={}
		longest_seq=[]
		filename=("".join([orthodir, '/', file]))
		orthogroup = os.path.splitext(file)[0]
		longdir=f"{separatedir}/{orthogroup}/{orthogroup}_longrep.fa"
		longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
		for record in SeqIO.parse(filename, "fasta"):
			id=record.id
			genome=id.split("_")[0]
			if genome not in protein_count:
				protein_count[genome]=0
			protein_count[genome]+=1
			if not longest or len(record) > len(longest):
				if genome not in longest.description:
					longest=record
		if longest:
			rename=SeqRecord(longest.seq, orthogroup, "", longest.id)

		counts=list(protein_count.values())
		equal=False
		if([counts[0]]*len(counts)==counts):
			equal=True
		if equal==False:
			blast_orthos.append(orthogroup)
			unique_counts=set(protein_count.values())
			blast_orthos_counts[orthogroup]=unique_counts
			os.makedirs(os.path.dirname(longdir),exist_ok=True)
			SeqIO.write(rename,longdir,'fasta')
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found ',str(len(blast_orthos)), ' orthologs to check.']))

else:
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	print("".join(['[', str(now), '] ','Found existing directory: ', separatedir, ' Continuing.']))

##Blast longest to all values
for orthogroup in os.listdir(separatedir):
	longfile=f"{separatedir}{orthogroup}/{orthogroup}_longrep.fa"
	seqfile=f"{orthodir}/{orthogroup}.fa"
	database=f"{separatedir}{orthogroup}/db/"
	os.makedirs(os.path.dirname(database), exist_ok=True)
	blastout=f"{separatedir}{orthogroup}/{orthogroup}_blastp.out"
	if not os.path.isfile(blastout):
		if os.path.isfile(longfile) and os.path.isfile(seqfile):
			shellfile=f"{separatedir}{orthogroup}/blastshell.sh"
			os.makedirs(os.path.dirname(shellfile),exist_ok=True)
			shell = open(shellfile, "w")
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
			shell.write("".join(['makeblastdb -in ', seqfile,' -dbtype prot -out ', database,
			orthogroup, '\n']))

			shell.write("".join(['blastp -query ', longfile, ' -db ', database, orthogroup,
			' -out ',blastout,
			' -outfmt 6 -evalue 1e-05', '\n']))

			shell.close()
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile, '\n']))
				os.system(execute)

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Completed BLASTs in: ', separatedir, ' Continuing.']))

##check if there are proteins aligning to the second portion of the longest sequence
similarity_check={}
for orthogroup in os.listdir(separatedir):
	aligned_genomes=[]
	blastout=f"{separatedir}{orthogroup}/{orthogroup}_blastp.out"
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	if os.path.isfile(blastout):
		potential_genomes=[]
		with open(blastout, 'r') as fin:
			lines=fin.readlines()
		for line in lines:
			sseqid=line.split("\t")[1]
			qstart=line.split("\t")[6]
			genome=sseqid.split("_")[0]
			if int(qstart)>50:
				potential_genomes.append(genome)

		for line in lines:
			sseqid=line.split("\t")[1]
			qstart=line.split("\t")[6]
			genome=sseqid.split("_")[0]
			if int(qstart)<=50 and genome in potential_genomes:
				aligned_genomes.append(genome)

	if len(aligned_genomes)>0:
		similarity_check[orthogroup]=aligned_genomes
		if not os.path.isdir(check_dir):
			os.makedirs(os.path.dirname(check_dir),exist_ok=True)

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Found ',str(len(similarity_check.keys())),
' orthogroups to check for overlaps and if neighboring each other. ']))

##Now find proteins that do not overlap AND are near eachother
no_overlap={}
ortho_prot_pairs={}
ortho_cuts={}
for orthogroup,genomes in similarity_check.items():
	genome_list=list(similarity_check[orthogroup])
	blastout=f"{separatedir}{orthogroup}/{orthogroup}_blastp.out"
	prot_pairs={}
	ortho_prots=[]
	if os.path.isfile(blastout):
		max_cuts=0
		for genome in genome_list:
			starts={}
			ends={}
			cuts=[]
			with open(blastout, 'r') as fin:
				lines=fin.readlines()
			for line in lines:
				sseqid=line.split("\t")[1]
				qstart=line.split("\t")[6]
				qend=line.split("\t")[7]
				genome_id=sseqid.split("_")[0]

				if genome_id==genome:
					starts[sseqid]=qstart
					ends[sseqid]=qend
			for protein1 in starts:
				start1=starts[protein1]
				end1=ends[protein1]
				if "-T" in protein1:
					prot1_num=int(((protein1.split("-T")[0]).split("_")[1]).lstrip("0"))
				else:
					prot1_num=int(((protein1.split("-t")[0]).split("_")[1]).lstrip("0"))
				for protein2 in starts:
					start2=starts[protein2]
					end2=ends[protein2]
					if "-T" in protein2:
						prot2_num=int(((protein2.split("-T")[0]).split("_")[1]).lstrip("0"))
					else:
						prot2_num=int(((protein2.split("-t")[0]).split("_")[1]).lstrip("0"))

					if prot2_num == (prot1_num +1) or prot2_num ==(prot1_num -1):
						if int(start2) > int(end1):
							if protein1 != protein2:
								prot_pairs[protein1]=protein2
								if protein1 not in cuts:
									cuts.append(protein1)
									ortho_prots.append(protein1)
								if protein2 not in cuts:
									cuts.append(protein2)
									ortho_prots.append(protein2)

			if len(cuts)>max_cuts:
				max_cuts=len(cuts)

	if max_cuts>0:
		ortho_cuts[orthogroup]=max_cuts
		no_overlap[orthogroup]=ortho_prots
	if len(prot_pairs)>0:
		ortho_prot_pairs[orthogroup]=prot_pairs

#maximum = max(ortho_cuts, key=ortho_cuts.get)
#print(maximum, ortho_cuts[maximum])

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Found ',str(len(no_overlap.keys())),' orthogroups with an average separation of ',
str(sum(ortho_cuts.values()) / len(ortho_cuts)),' to create queries for. ']))

#Now put representatives in a file for blast
blast_pairs={}
for orthogroup in no_overlap:
	protein_list=list(no_overlap[orthogroup])
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	blast_seqs=f"{check_dir}{orthogroup}_queries.fa"
	protein_seqs=[]
	if not os.path.isfile(blast_seqs):
		if os.path.isdir(check_dir):
			filename=f"{orthodir}/{orthogroup}.fa"
			if os.path.isfile(filename):
				for record in SeqIO.parse(filename, "fasta"):
					id=record.id
					if id in protein_list:
						protein_seqs.append(record)
		SeqIO.write(protein_seqs,blast_seqs,'fasta')

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Found ',str(len(no_overlap.keys())),
' orthogroups with no-overlap representatives to check for similarities.']))

#blast overlaps
for orthogroup in no_overlap:
	database=f"{separatedir}{orthogroup}/db/"
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	blast_seqs=f"{check_dir}{orthogroup}_queries.fa"
	blastout=f"{check_dir}{orthogroup}_blastp.out"
	if not os.path.isfile(blastout):
		if os.path.isfile(blast_seqs):
			shellfile=f"{check_dir}blastshell.sh"
			os.makedirs(os.path.dirname(shellfile),exist_ok=True)
			shell=open(shellfile,"w")
			shell = open(shellfile, "w")
			shell.write("""#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=""")
			shell.write(f"{check_dir}")
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
			shell.write("".join(['cd ', check_dir, '\n', '\n']))

			#run blast
			shell.write("".join(['blastp -query ', blast_seqs, ' -db ', database, orthogroup,
			' -out ',blastout,' -outfmt 6', '\n']))

			shell.close()
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile, '\n']))
				os.system(execute)

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Completed BLASTs in: ', str(len(no_overlap.keys())),
' represenatives. Continuing.']))

#start grouping
ortholog_groups={}
for orthogroup in no_overlap:
	groups={}
	group_count=0
	unique_groups=[]
	protein_list=list(no_overlap[orthogroup])
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	blastout=f"{check_dir}{orthogroup}_blastp.out"
	if os.path.isfile(blastout):
		for protein1 in protein_list:
			if protein1 not in unique_groups:
				unique_groups.append(protein1)
				protein1_group=[]
				protein1_group.append(protein1)
				with open(blastout,'r') as fin:
					lines=fin.readlines()
				for line in lines:
					qseqid =line.split("\t")[0]
					if protein1 == qseqid:
						protein2=line.split("\t")[1]
						evalue=line.split("\t")[10]
						if protein2 in protein_list:
							if protein2 not in unique_groups:
								if Decimal(evalue) < Decimal('1.0e-50'):
									if protein2 not in protein1_group:
										protein1_group.append(protein2)
										unique_groups.append(protein2)
				group_count+=1
				groups[f"group_{group_count}"]=protein1_group
	ortholog_groups[orthogroup]=groups

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Found ',str(len(ortholog_groups.keys())),
' orthogroups to split into an average of ', str(sum(ortho_cuts.values()) / len(ortho_cuts)),' groups.']))

#Write longest representative for each group into a file
group_reps={}
for orthogroup in ortholog_groups:
	groups=ortholog_groups[orthogroup]
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	seqs_file=f"{check_dir}{orthogroup}_queries.fa"
	longdir=f"{separatedir}{orthogroup}/split/{orthogroup}_longrep.fa"
	os.makedirs(os.path.dirname(longdir), exist_ok=True)
	longest_group_reps=[]
	group_ids=[]
	if os.path.isfile(seqs_file):
		for group in groups:
			longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			reps=list(groups[group])
			for record in SeqIO.parse(seqs_file, "fasta"):
				if record.id in reps:
					if not longest or len(record) > len(longest):
						if record.id not in group_ids:
							longest=record
			if longest:
				group_ids.append(longest.id)
				rename=SeqRecord(longest.seq,group,"",longest.id)
				longest_group_reps.append(rename)

	SeqIO.write(longest_group_reps,longdir,'fasta')
	group_reps[orthogroup]=group_ids

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Created ',str(len(group_reps.keys())),
' groups for combined ortholog representatives.' ]))

#now blast all groups
for orthogroup in group_reps:
	check_dir=f"{separatedir}{orthogroup}/similarity_check/"
	longfile=f"{separatedir}{orthogroup}/split/{orthogroup}_longrep.fa"
	blastout=f"{separatedir}{orthogroup}/split/{orthogroup}_group.out"
	subject_seqs=f"{orthodir}/{orthogroup}.fa"
	database=f"{separatedir}{orthogroup}/db/"
	if not os.path.isfile(blastout):
		if os.path.isfile(longfile) and os.path.isfile(subject_seqs):
			shellfile=f"{separatedir}{orthogroup}/split/blastshell.sh"
			os.makedirs(os.path.dirname(shellfile),exist_ok=True)
			shell=open(shellfile,"w")
			shell = open(shellfile, "w")
			shell.write("""#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=""")
			shell.write(f"{separatedir}{orthogroup}/split/")
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
			shell.write("".join(['makeblastdb -in ', longfile,' -dbtype prot -out ', database,
			orthogroup, '_groups','\n']))

			#THIS CHANGES THE NORMAL OUTPUT
			shell.write("".join(['blastp -query ', subject_seqs, ' -db ', database, orthogroup,
			'_groups -out ',blastout,
			' -outfmt "6 qseqid sseqid qstart qend qcovs evalue bitscore" -evalue 1e-05', '\n']))

			shell.close()
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile, '\n']))
				os.system(execute)

now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print("".join(['[', str(now), '] ','Completed BLASTs. Continuing.']))

##Pick out proteins to split
split_proteins={}
CNAG_split=[]
orthogroup_split=[]
for orthogroup in group_reps:
	blastout=f"{separatedir}{orthogroup}/split/{orthogroup}_group.out"
	split_file=f"{separatedir}{orthogroup}/{orthogroup}_splits.fa"
	full_seqs=f"{orthodir}/{orthogroup}.fa"
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
			evalue=line.split("\t")[5]
			bitscore=(line.split("\t")[6]).split("\n")[0]
			if int(qcovs) < 90:
				if Decimal(evalue) < Decimal('1.0e-100'):
					if qseqid not in protein_count:
						protein_count[qseqid]=0
					protein_count[qseqid]+=1
		for protein in protein_count:
			split_sites=[]
			if protein_count[protein]>1:
				groups=[]
				for line in lines:
					if protein in line:
						sseqid=line.split("\t")[1]
						qstart=line.split("\t")[2]
						qend=line.split("\t")[3]
						groups.append(sseqid)
				for group in groups:
					for line in lines:
						if group and protein in line:
							qstart=line.split("\t")[2]
							qend=line.split("\t")[3]
							evalue=line.split("\t")[5]
							if Decimal(evalue) < Decimal('1.0e-100'):
								if int(qend)-int(qstart)>=35:
									split_site=f"{qstart}:{qend}"
									if split_site not in split_sites:
										split_sites.append(split_site)
										if "CNAG" in protein:
											if orthogroup not in CNAG_split:
												CNAG_split.append(orthogroup)
			if len(split_sites)>0:
				split_proteins[protein]=split_sites
				if orthogroup not in orthogroup_split:
					orthogroup_split.append(orthogroup)
					print(orthogroup)
print(str(len(CNAG_split)))
print(str(len(orthogroup_split)))

#write new split fastas into a new file
all_split_count=0
numbers = [1,2,3,4,5]
alphabet = string.ascii_uppercase
mapping = {number: alphabet[i] for i, number in enumerate(numbers)}
new_prots_dir=f"{userdir}/new_protein_files/"
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
print("".join(['[', str(now), '] ','Split ', str(all_split_count), ' proteins. New protein file located in: ', new_prots_dir]))
