##Last updated: 6/12/26

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
parser = argparse.ArgumentParser(description='Submodule of miniprot.py that renames proteins GACRC Slurm',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
args = parser.parse_args()

##file definitions
dir = os.getcwd()
directory_header=args.dir[0:10].replace("_","-")
current_orthos=f"{dir}/{directory_header}.csv"
sorted_gff=f"{dir}/{args.dir}_minisorted.gff3"
old_prots=f"{dir}/{args.dir}_miniprot.fa"
final_gff=f"{dir}/{args.dir}.gff3"
new_prots=f"{dir}/{args.dir}_newprots.fa"

##get last annot number
if os.path.isfile(sorted_gff):
	last_prot=[]
	with open(sorted_gff,'r') as fin:
		for line in fin:
			if "#" not in line:
				att=line.split("\t")[8]
				if ";" in att:
					id=att.split(";")[0]
				else:
					id=att
				if "OG" not in id:
					prot_num=(id.split("=")[1]).split("_")[1]
					if "-T" in prot_num:
						only_number=prot_num.split("-T")[0]
					else:
						only_number=prot_num.split("\n")[0]
					if only_number not in last_prot:
						last_prot.append(int(only_number))
	prot_num=int(max(last_prot))

##Go through and replace orthogroups with new protein ids
og_to_prot={}
last_protnum=prot_num
if os.path.isfile(sorted_gff):
	seen=set()
	with open(final_gff,'w') as fout:
		with open(sorted_gff,'r') as fin:
			for line in fin:
				if "OG000" in line:
					att=line.split("\t")[8]
					id=att.split(";")[0]
					id_num=id.split("=")[1]
					ortho=id.split("M")[0]
					if id_num not in og_to_prot:
						last_protnum+=1
						new_id=f"{directory_header}_{last_protnum:06d}-T1"
						og_to_prot[id_num]=new_id
						new_line=line.replace(id_num,new_id)
						if new_line not in seen:
							fout.write(new_line)
							seen.add(new_line)
					if id_num in og_to_prot:
						new_id=og_to_prot[id_num]
						new_line=line.replace(id_num,new_id)
						if new_line not in seen:
							fout.write(new_line)
							seen.add(new_line)
				else:
					if "###" not in line:
						if line not in seen:
							fout.write(line)
							seen.add(line)

##Now write final protein names to file
if os.path.isfile(old_prots):
	new_names=[]
	unique=[]
	for record in SeqIO.parse(old_prots,"fasta"):
		miniid=(record.id).split("_")[1]
		if miniid in og_to_prot:
			protid=og_to_prot[miniid]
			seq_fa=SeqRecord(record.seq,protid,"",miniid)
			if protid not in unique:
				new_names.append(seq_fa)
				unique.append(protid)
	sorted_new=sorted(new_names,key=lambda r:r.id)
	SeqIO.write(sorted_new,new_prots,"fasta")
