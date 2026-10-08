##Last updated: 6/1/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
parser = argparse.ArgumentParser(description='submodule of miniprot.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
parser.add_argument('ortho_dir')
args = parser.parse_args()

##file definitions
dir = os.getcwd()

##get all missing orthologs
directory_header=args.dir[0:10].replace("_","-")
current_orthos=f"{dir}/{directory_header}.csv"
missing_orthos=f"{dir}/{args.dir}_missing_orthos.fa"
if not os.path.isfile(current_orthos):
	current_orthologs={}
	longest_seqs=[]
	for file in os.listdir(args.ortho_dir):
		proteins_found=[]
		orthogroup=os.path.splitext(file)[0]
		filename=f"{args.ortho_dir}/{file}"
		for record in SeqIO.parse(filename,"fasta"):
			genome = (record.id).split("_")[0]
			if directory_header == genome:
				proteins_found.append(record.description)
				if record.description not in current_orthologs:
					current_orthologs[record.description]=orthogroup

		#write sequences of missing orthologs
		if len(proteins_found)==0:
			longest_non_CNAG=SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			longest_CNAG=SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			for record in SeqIO.parse(filename,"fasta"):
				genome=(record.id).split("_")[0]
				if "CNAG" in genome:
					if not longest_CNAG or len(record)> len(longest_CNAG):
						if genome not in longest_CNAG.description:
							longest_CNAG=record
				else:
					if not longest_non_CNAG or len(record)>len(longest_non_CNAG):
						if genome not in longest_non_CNAG.description:
							longest_non_CNAG=record
			if "CNAG" in longest_CNAG.description:
				rename_CNAG=SeqRecord(longest_CNAG.seq,orthogroup,"",longest_CNAG.id)
				longest_seqs.append(rename_CNAG)
			else:
				rename_non_CNAG=SeqRecord(longest_non_CNAG.seq,orthogroup,"",longest_non_CNAG.id)
				longest_seqs.append(rename_non_CNAG)

	SeqIO.write(longest_seqs,missing_orthos,'fasta')

	#write to csv
	with open(current_orthos,'w') as outfile:
		for prot,ortho in current_orthologs.items():
			outfile.write("".join([ortho,',',prot,'\n']))

	if len(current_orthologs)==0:
		print("".join(['Unable to find any orthologs for header: ', directory_header,
				'Please check this is the correct protein heading for genome: ',
					args.dir]))
else:
	print("".join(["Found current orthologs under ",current_orthos," Skipping get-orthos.py"]))
