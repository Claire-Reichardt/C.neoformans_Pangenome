##Last updated:	3/17/25

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
parser = argparse.ArgumentParser(description='Script for removing bad Orthologs based on read depth comparison on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)', required=True)
required.add_argument('-b',dest='buscodir', metavar='BUSCO_DIRECTORY', help='Directory containing protein fastas for BUSCO single copy orthologs from busco_map.py', required=True)
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
readdir = f"{userdir}/readmapping"
removefile = f"{userdir}/bad_orthologs.tsv"
addfile = f"{userdir}/missing_orthologs.tsv"

#Check directories exist
if os.path.isdir(readdir) and os.access(readdir, os.R_OK):
	readdir = (os.path.abspath(readdir))
else:
	sys.exit(''.join(["ERROR: readmapping directory:", readdir, " is missing or not readable!!"]))

GENOMESDIR =f"{args.orthodir}/all_orthos"
if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")

BUSCODIR = args.buscodir
if os.path.isdir(BUSCODIR) and os.access(BUSCODIR, os.R_OK):
	buscodir = (os.path.abspath(BUSCODIR))
else:
	sys.exit("ERROR: BUSCO sequence directory is missing or not readable! Please make sure to execute busco_map.py")

#Find orthologs with multiple maps (duplicates)
ortho_mapq={}
ortho_counts={}
for directory in os.listdir(readdir):
	strain=directory
	straindir=f"{userdir}/readmapping/{directory}"
	if os.path.isdir(straindir):
		covfile=f"{straindir}/{strain}_coverage.tsv"
		if os.path.isfile(covfile):
			with open(covfile,'r') as fin:
				lines=fin.readlines()
			for line in lines:
				if "#rname" not in line:
					ortho=line.split("\t")[0]
					average=line.split("\t")[6]
					coverage=line.split("\t")[5]
					mapq=(line.split("\t")[8]).split("\n")[0]
					if ortho not in ortho_mapq:
						ortho_mapq[ortho]=0
						ortho_counts[ortho]=0
					ortho_mapq[ortho]+=int(float(mapq))
					ortho_counts[ortho]+=1
multi_mappers=[]
for ortholog in ortho_mapq:
	total_mapq=ortho_mapq[ortholog]
	count=ortho_counts[ortholog]
	average_mapq=total_mapq/count
	if average_mapq < 10:
		if ortholog not in multi_mappers:
			multi_mappers.append(ortholog)
			print(average_mapq)
			#print(ortholog)
			#print(count)
print(int(len(multi_mappers)))
print(multi_mappers)
#Make a key of orthologs with associated proteins
prot_ortho_dict={}
for file in os.listdir(genomesdir):
	filename=os.path.join(genomesdir,file)
	with open(filename) as fin:
		lines=fin.readlines()
	for line in lines:
		if "\t" in line:
			ortholog=line.split("\t")[0]
			protein=(line.split("\n")[0]).split("\t")[1]
			prot_ortho_dict[protein]=ortholog

#Get list of single copy orthologs (SCOs) from BUSCO
SCO_dict={}
multi_SCO={}
for file in os.listdir(buscodir):
	SCO_ortho_counter={}
	SCO_ortho_list=[]
	strain=file.split("_SCO.fa")[0]
	filename=os.path.join(buscodir,file)
	for record in SeqIO.parse(filename,"fasta"):
		id=record.id
		if id in prot_ortho_dict:
			ortholog=prot_ortho_dict[id]
			if ortholog not in SCO_ortho_counter:
				SCO_ortho_counter[ortholog]=0
			SCO_ortho_counter[ortholog]+=1
	for ortholog in SCO_ortho_counter:
		if SCO_ortho_counter[ortholog]==1:
			SCO_ortho_list.append(ortholog)
	SCO_dict[strain]=SCO_ortho_list
prot_ortho_dict.clear()

#Go through every strain directory
dict_check={}
dict_remove={}
for directory in os.listdir(readdir):
	check_orthos=[]
	remove_orthos=[]
	strain=directory
	straindir=f"{userdir}/readmapping/{directory}"
	if os.path.isdir(straindir):
		covfile=f"{straindir}/{strain}_coverage.tsv"
		if os.path.isfile(covfile):
			if strain in SCO_dict:
				SCOs=SCO_dict[strain]
				dict_medorthos={}
				dict_coverage={}
				dict_qual={}
				averages=[]
				with open(covfile,'r') as fin:
					lines=fin.readlines()
				for line in lines:
					if "#rname" not in line:
						ortho=line.split("\t")[0]
						average=line.split("\t")[6]
						coverage=line.split("\t")[5]
						mapq=line.split("\t")[8]
						if ortho in SCOs:
							averages.append(int(float(average)))
						dict_coverage[ortho]=int(float(coverage))
						dict_qual[ortho]=int(float(mapq))
						if int(float(average)) !=0:
							dict_medorthos[ortho]=int(float(average))
						if int(float(average))==0:
							dict_medorthos[ortho]=0.001
				med_reads=statistics.median(averages)
				dict_norm_reads={key: value/(med_reads) for key, value in dict_medorthos.items()}
				dict_log_norm={key: np.log2(value) for key, value in dict_norm_reads.items()}
				for key,value in dict_log_norm.items():
					cover=dict_coverage[key]
					mapq=dict_qual[key]
					if value >= 0:
						#if mapq >=45:
						if cover >= 50:
							check_orthos.append(key)
					if value <= -5:
						#if mapq <=45:
						#	if cover <= 50:
						remove_orthos.append(key)

		#print if there is a problem finding the BUSCO SCOs
			else:
				print("".join(['BUSCO SCO file not found for: ', strain, ' in: ', buscodir]))
		#print if there is a problem finding the average depth file
		else:
			print("".join(['Read depth analysis did not run correctly for: ', straindir]))

	#put orthos to check in a dictonary
	dict_check[strain]=check_orthos
	dict_remove[strain]=remove_orthos

SCO_dict.clear()

#Compare present orthos to big hits from read mapping
fake_orthos=[]
added_count=0
for key,value in dict_check.items():
	to_check=list(dict_check[key])
	strain_file=f"{genomesdir}/{key}.tsv"
	added_orthos=[]
	if os.path.isfile(strain_file):
		present_orthos=[]
		with open(strain_file,'r') as fin:
			lines=fin.readlines()
		for line in lines:
			ortho=line[0:9]
			present_orthos.append(ortho)
		for missing_ortho in to_check:
			if missing_ortho not in present_orthos:
				if missing_ortho not in added_orthos:
					with open(strain_file, 'a') as file:
						file.write("".join([missing_ortho, '\n']))
						added_count+=1
						added_orthos.append(missing_ortho)
					if missing_ortho not in fake_orthos:
						fake_orthos.append(missing_ortho)

	else:
		print("".join(['Unable to find ', strain_file, ' in ', genomesdir]))

#Remove orthologs with poor presence from each strain's ortholog amount
removed_orthos=[]
removed_count=0
for strain,ortho in dict_remove.items():
	to_remove=list(dict_remove[strain])
	strain_file=f"{genomesdir}/{strain}.tsv"
	if os.path.isfile(strain_file):
		print(strain_file)
		with open(strain_file, 'r') as input:
			with open("temp.txt",'a') as output:
				for line in input:
					ortho_check=line[0:9]
					if ortho_check not in to_remove:
						if ortho_check in line.strip("\n"):
							output.write(line)
					if ortho_check in to_remove:
						if ortho_check in line.strip("\n"):
							removed_count+=1
							if ortho_check not in removed_orthos:
								removed_orthos.append(ortho_check)

		os.replace("temp.txt", strain_file)

#Print how many fake orthos were found
print("".join(['Found ', str(added_count),' inconsistant orthologs that were appended']))
print("".join(['Found ', str(removed_count), ' inconsistant orthologs that were removed']))

#Write orthologs to a file
n_fake_orthos= ["{}\n".format(i) for i in sorted(fake_orthos)]
with open(addfile, 'w') as fp:
	fp.writelines(n_fake_orthos)

n_removed_orthos= ["{}\n".format(i) for i in sorted(removed_orthos)]
with open(removefile, 'w') as fp:
        fp.writelines(n_removed_orthos)

