##Last updated: 8/3/26

##import system
import sys
import subprocess
import os
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
from pathlib import Path
import shutil
parser = argparse.ArgumentParser(description='Script to complete miniprot.py data',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY_DIR", dest="asbdir", help="Directory containing list of long read assembilies", required=True)
required.add_argument("-m", metavar="MINI_DIR", dest="minidir", help="miniprot directory created by miniprot.py", required=True)
required.add_argument('-c', metavar='ORTHOLOG_CLASSIFICATION', dest="classify",help='classify.csv from new_classify.py output')

args = parser.parse_args()

##Check files exist
if os.path.isdir(args.minidir) and os.access(args.minidir, os.R_OK):
	minidir=(os.path.abspath(args.minidir))
else:
	sys.exit("ERROR: miniprot directory is missing or not readable!!")

CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.csv is missing or not readable!!")

if os.path.isdir(args.asbdir) and os.access(args.asbdir, os.R_OK):
	asbdir=(os.path.abspath(args.asbdir))
else:
	sys.exit("ERROR: Assembly directory is missing or not readable!!")

##file definitions
userdir = os.getcwd()
sub_out=f"{userdir}/subtelo_stats.csv"

#Get the classification for each ortholog
class_dict={}
if os.path.isfile(classfile):
	with open(classfile,'r') as fin:
		lines=fin.readlines()
	for line in lines:
		no_n=line.split("\n")[0]
		ortho=line[0:9]
		gene_class=no_n.split(",")[1]
		class_dict[ortho]=gene_class

##Go through each directory
per_subtelo_acc={}
for file in os.listdir(asbdir):
	filename=("".join([asbdir,'/',file]))
	genome_name=os.path.splitext(file)[0]
	header=genome_name[0:10].replace('_','-')
	genome_dir=f"{minidir}/{genome_name}"
	if os.path.isdir(genome_dir):
		prot_dict={}
		scaffold_size={}
		acc_orthos=[]
		acc_subtelo=[]
		genome_size=0
		aln_prots=f"{genome_dir}/{genome_name}_aln.gff"
		current_orthos=f"{genome_dir}/{header}.csv"
		mapped_ortho_file=f"{genome_dir}/{genome_name}_mapped_orthos.csv"

		#get all current proteins
		with open(current_orthos,'r') as fin:
			for line in fin:
				ortho=line.split(",")[0]
				if "M" in ortho:
					ortho=(line.split(",")[0]).split("M")[0]
				protein=(line.split(",")[1]).split("\n")[0]
				if " " in protein:
					protein=line.split(",")[1].split(" ")[0]
				prot_dict[protein]=ortho

		#get all chromosome lengths
		for record in SeqIO.parse(filename,"fasta"):
			scaffold_size[record.id]=len(record)
			genome_size+=len(record)

		#only look at assemblies with 14 chromosomes
		if len(scaffold_size)==14:
			with open(aln_prots,'r') as fin:
				for line in fin:
					if "##" not in line:
						feature=line.split("\t")[2]
						if feature=="mRNA":
							att=(line.split("\t")[8]).split("\n")[0]
							protein=(att.split("Target=")[1]).split(" ")[0]
							if "OG" not in protein and protein in prot_dict:
								ortho=prot_dict[protein]
								oclass=class_dict[ortho]
								if oclass=="accessory":
									chrom=line.split("\t")[0]
									start=int(line.split("\t")[3])
									stop=int(line.split("\t")[4])
									ident=float((att.split(";")[2]).split("Identity=")[1])
									if ident>0.85:
										chr_length=scaffold_size[chrom]
										chr_end=int(chr_length)-40000
										if start < 40000 or stop > chr_end:
											if protein not in acc_subtelo:
												acc_subtelo.append(protein)
										else:
											if protein not in acc_orthos:
												acc_orthos.append(protein)
			#calculate percents
			total_acc=len(acc_subtelo)+len(acc_orthos)
			percent=int(len(acc_subtelo))/total_acc
			with open(sub_out,"a") as f:
				f.write(''.join([genome_name,",",str(genome_size),",",str(len(acc_subtelo)),",",
					str(total_acc),",",str(percent),"\n"]))
			print(''.join(['Found ', str(percent), ' accessory genes in the subtelomeres of ', genome_name]))
		else:
			print("".join([genome_name," is not 14 scaffolds so it will be ignored",]))

	else:
		sys.exit("".join(['ERROR : Could not find miniprot directory for ',genome_name,' in ',minidir]))


