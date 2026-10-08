##Last updated: 6/16/26

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
parser = argparse.ArgumentParser(description='Submodule of miniprot.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
parser.add_argument('fun_gff')
args = parser.parse_args()

##file definitions
dir = os.getcwd()
longdir=os.path.dirname(os.path.dirname(dir))
longfile=f"{longdir}/longest_missing_orthos.fa"

##Search completed
gff_file=f"{dir}/{args.dir}_aln.gff"
hit_gff=f"{dir}/{args.dir}_hits.gff"
directory_header=args.dir[0:10].replace('_','-')
current_orthos=f"{dir}/{directory_header}.csv"
fun_gff=f"{args.fun_gff}/cryptococcus_{args.dir}.gff3"
new_gff=f"{dir}/{args.dir}_miniprot.gff3"
overlap_file=f"{dir}/{args.dir}_overlapping-orthologs.csv"
new_prots=f"{dir}/{args.dir}_miniprot.fa"
##Make a dictionary of proteins and their orthogroups
prot_dict={}
parent_names=[]
if os.path.isfile(current_orthos):
	with open(current_orthos,'r') as fin:
		for line in fin:
			protein=(line.split(",")[1]).split("\n")[0]
			parent=protein.split("-T1")[0]
			orthogroup=line.split(",")[0]
			if protein not in prot_dict:
				prot_dict[protein]=orthogroup
			if parent not in parent_names:
				parent_names.append(parent)

##Search for good hits
if os.path.isfile(current_orthos) and os.path.isfile(gff_file):
	best_hits=[]
	with open(gff_file,'r') as fin:
		for line in fin:
			if "mRNA" in line:
				att=line.split("\t")[8]
				id=att.split(";")[0]
				id_num=id.split("ID=")[1]
				ortho=id.split("M")[0]
				perid=att.split(";")[2]
				perid_only=float(perid.split("Identity=")[1])
				if perid_only > 0.75:
					best_hits.append(id_num)

	with open(hit_gff,"w") as f:
		for id_num in best_hits:
			ortholog=id_num.split("M")[0]
			with open(gff_file,'r') as fin:
				for line in fin:
					if ortholog in line:
						f.write(line)

##Now see if CNAG protein hits overlap with previous gff
CNAG_orthos=[]
for record in SeqIO.parse(longfile,"fasta"):
	if "CNAG" in record.description:
		if record.id not in CNAG_orthos:
			CNAG_orthos.append(record.id)

def is_overlapping(x1,x2,y1,y2):
	return  max(x1,y1) < min(x2,y2)

new_orthoids=[]

#go through and get all current gene start and stops
gene_locs={}

with open(fun_gff,'r') as fin:
	for line in fin:
		if "mRNA" in line:
			scaffold=line.split("\t")[0]
			f_start=line.split("\t")[3]
			f_end=line.split("\t")[4]
			id_num=((line.split("\t")[8]).split(";")[0]).split("ID=")[1]
			range=f"{id_num}:{f_start}-{f_end}"
			if scaffold in gene_locs:
				gene_locs[scaffold].append(range)
			else:
				gene_locs[scaffold]=[range]

fin.close()

#map all cnag hits and update gene locations
overlapping_orthos=[]
with open(hit_gff,'r') as fin:
	for line in fin:
		if "mRNA" in line:
			scaffold=line.split("\t")[0]
			m_start=int(line.split("\t")[3])
			m_end=int(line.split("\t")[4])
			att=line.split("\t")[8]
			id=att.split(";")[0]
			id_num=id.split("ID=")[1]
			ortho=id_num.split("M")[0]

			if ortho in CNAG_orthos:
				overlap_list=[]
				overlap_orthos=[]

				if scaffold in gene_locs:
					for location in gene_locs[scaffold]:
						id2=location.split(":")[0]
						f_start=int((location.split(":")[1]).split("-")[0])
						f_end=int((location.split(":")[1]).split("-")[1])
						overlap=is_overlapping(m_start,m_end,f_start,f_end)
						if overlap==True:
							if id2 in prot_dict:
								overortho=prot_dict[id2]
								if overortho != ortho:
									overlap_orthos.append(overortho)
							if overlap not in overlap_list:
								overlap_list.append(overlap)
								#print("".join([str(overlap)," ",str(m_start), " ",str(m_end)," ",str(f_start)," ",str(f_end)]))
				else:
					range=f"{id_num}:{m_start}-{m_end}"
					gene_locs[scaffold]=[range]

				if True not in overlap_list:
					if id_num not in new_orthoids:
						new_orthoids.append(id_num)
						range=f"{id_num}:{m_start}-{m_end}"
						if scaffold in gene_locs:
							gene_locs[scaffold].append(range)
						else:
							gene_locs[scaffold]=[range]

				if len(overlap_orthos)>0:
					for item in overlap_orthos:
						ortho_ids=f"{ortho}:{item}"
						if ortho_ids not in overlapping_orthos:
							overlapping_orthos.append(ortho_ids)
fin.close()

#now finally map non-cnag orthologs and append locations
with open(hit_gff,'r') as fin:
	for line in fin:
		if "mRNA" in line:
			scaffold=line.split("\t")[0]
			m_start=int(line.split("\t")[3])
			m_end=int(line.split("\t")[4])
			att=line.split("\t")[8]
			id=att.split(";")[0]
			id_num=id.split("ID=")[1]
			ortho=id_num.split("M")[0]

			if ortho not in CNAG_orthos:
				overlap_list=[]
				overlap_orthos=[]

				if scaffold in gene_locs.keys():
					for location in gene_locs[scaffold]:
						id2=location.split(":")[0]
						f_start=int((location.split(":")[1]).split("-")[0])
						f_end=int((location.split(":")[1]).split("-")[1])
						overlap=is_overlapping(m_start,m_end,f_start,f_end)
						if overlap==True:
							if id2 in prot_dict:
								overortho=prot_dict[id2]
								if overortho != ortho:
									overlap_orthos.append(overortho)
							if overlap not in overlap_list:
								overlap_list.append(overlap)
								#print("".join([str(overlap)," ",str(m_start), " ",str(m_end)," ",str(f_start)," ",str(f_end)]))
				else:
					range=f"{id_num}:{m_start}-{m_end}"
					gene_locs[scaffold]=[range]

				if True not in overlap_list:
					if id_num not in new_orthoids:
						new_orthoids.append(id_num)
						range=f"{id_num}:{m_start}-{m_end}"
						if scaffold in gene_locs:
							gene_locs[scaffold].append(range)
						else:
							gene_locs[scaffold]=[range]

				if len(overlap_orthos)>0:
					for item in overlap_orthos:
						ortho_ids=f"{ortho}:{item}"
						if ortho_ids not in overlapping_orthos:
							overlapping_orthos.append(ortho_ids)

fin.close()

#write a file of orthologs that overlap with cnags for possible recombining
missing_orthos=[]
for record in SeqIO.parse(longfile,"fasta"):
	if record.id not in missing_orthos:
		missing_orthos.append(record.id)
with open(overlap_file,"w") as f:
	for item in overlapping_orthos:
		ortho1=item.split(":")[0]
		ortho2=item.split(":")[1]
		if ortho2 in missing_orthos:
			f.write("".join([ortho1,",",ortho2,"\n"]))
f.close()

##write good orthologs to gff and protein files
new_seqs=[]
out_gff=open(new_gff,"w")
with open(fun_gff,"r") as fin:
	for line in fin:
		out_gff.write(line)

fin.close()

#only write proteins >50aa for non-CNAGs
good_orthologs=[]
for orthoid in new_orthoids:
	ortholog=orthoid.split("M")[0]
	with open(gff_file,"r") as fin:
		prevLine=""
		for line in fin:
			if "mRNA" in line:
				if orthoid in line:
					sta=prevLine
					seq=Seq((sta.split("##STA\t")[1]).split("\n")[0])
					if ortholog in CNAG_orthos:
						good_orthologs.append(orthoid)
					else:
						if len(seq)>50:
							good_orthologs.append(orthoid)
			prevLine=line
fin.close()

for orthoid in good_orthologs:
	ortholog=orthoid.split("M")[0]
	with open(hit_gff,"r") as fin:
		for line in fin:
			if orthoid in line:
				out_gff.write(line)

	with open(gff_file,"r") as fin:
		prevLine=""
		for line in fin:
			if "mRNA" in line:
				if orthoid in line:
					sta=prevLine
					seq=Seq((sta.split("##STA\t")[1]).split("\n")[0])
					header=f"{directory_header}_{orthoid}"
					seq_fa=SeqRecord(seq,header,"Hit from miniprot",orthoid)
					new_seqs.append(seq_fa)
			prevLine=line

SeqIO.write(new_seqs,new_prots,"fasta")

print("".join(['Found ',str(len(good_orthologs)),' new orthologs for ',dir]))
out_gff.close()
fin.close()
