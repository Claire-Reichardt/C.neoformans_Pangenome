##Last updated: 7/7/26

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
args = parser.parse_args()

##file definitions
dir = os.getcwd()

##Search completed
gff_file=f"{dir}/{args.dir}_aln.gff"
missing_file=f"{dir}/{args.dir}_missing.gff"
hit_gff=f"{dir}/{args.dir}_hits.gff"
directory_header=args.dir[0:10].replace('_','-')
current_orthos=f"{dir}/{directory_header}.csv"
missing_prots=f"{dir}/{args.dir}_missing.fa"
missing_orthos=f"{dir}/{args.dir}_missing_orthos.csv"
mapped_orthos=f"{dir}/{args.dir}_mapped_orthos.csv"

##Make a dictionary of proteins and their orthogroups
prot_dict={}
if os.path.isfile(current_orthos):
	with open(current_orthos,'r') as fin:
		for line in fin:
			protein=(line.split(",")[1]).split("\n")[0]
			if " " in protein:
				protein=((line.split(",")[1]).split("\n")[0]).split(" ")[0]
			print(protein)
			orthogroup=line.split(",")[0]
			if protein not in prot_dict:
				prot_dict[protein]=orthogroup

##Search for good hits
if os.path.isfile(current_orthos) and os.path.isfile(missing_file):
	best_hits=[]
	with open(missing_file,'r') as fin:
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
			with open(missing_file,'r') as fin:
				for line in fin:
					if ortholog in line:
						f.write(line)
	print("".join(['Found ', str(len(best_hits)), ' to check for overlaps']))

##define overlapping hits
def is_overlapping(x1,x2,y1,y2):
	return  max(x1,y1) < min(x2,y2)

##map all hits and gene locations
gene_locs={}
found_orthos=[]
with open(gff_file,'r') as fin:
	for line in fin:
		if "mRNA" in line:
			scaffold=line.split("\t")[0]
			f_start=line.split("\t")[3]
			f_end=line.split("\t")[4]
			id_num=((line.split("\t")[8]).split(";")[0]).split("ID=")[1]
			prot=id_num.rsplit("M",1)[0]
			range=f"{id_num}:{f_start}-{f_end}"
			if scaffold in gene_locs:
				gene_locs[scaffold].append(range)
			else:
				gene_locs[scaffold]=[range]
			if prot in prot_dict:
				ortholog=prot_dict[prot]
				if ortholog not in found_orthos:
					found_orthos.append(ortholog)
fin.close()

found_orthos.sort()
out_found=open(mapped_orthos,"w")
for ortholog in found_orthos:
	out_found.write("".join([ortholog,"\n"]))

##see if hits are overlapping
overlapping_orthos=[]
new_orthoids=[]
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
			overlap_list=[]
			overlap_orthos=[]

			if scaffold in gene_locs:
				for location in gene_locs[scaffold]:
					id2=(location.split(":")[0]).split("M")[0]
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

##Only look for proteins >50aa
good_orthologs=[]
for orthoid in new_orthoids:
	ortholog=orthoid.split("M")[0]
	with open(missing_file,'r') as fin:
		prevLine=""
		for line in fin:
			if "mRNA" in line:
				if orthoid in line:
					sta=prevLine
					seq=Seq((sta.split("##STA\t")[1]).split("\n")[0])
					if len(seq)>50:
						if orthoid not in good_orthologs:
							good_orthologs.append(orthoid)
							#print("".join([orthoid,":",str(len(seq))]))
			prevLine=line
fin.close()

##write missing orthologs
final_out=[]
good_orthologs.sort()
out_orthos=open(missing_orthos,"w")
for orthoid in good_orthologs:
	ortholog=orthoid.split("M")[0]
	if ortholog not in final_out:
		out_orthos.write("".join([ortholog,"\n"]))
		final_out.append(ortholog)
