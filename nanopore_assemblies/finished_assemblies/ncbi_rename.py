##Last updated: 6/9/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Code to edit fasta headers for ncbi submission',
epilog='Written by Claire Reichardt')

parser.add_argument('assembly_file')
parser.add_argument('isolate_name')
args=parser.parse_args()

dir=os.getcwd()
ncbi_dir=f"{dir}/ncbi_final/"
os.makedirs(os.path.dirname(ncbi_dir),exist_ok=True)
new_file=f"{ncbi_dir}{args.isolate_name}.fsa"

if os.path.isfile(args.assembly_file):
	with open(new_file,'a') as out_f:
		asb_file=os.path.abspath(args.assembly_file)
		chrom=0
		for record in SeqIO.parse(asb_file,"fasta"):
			chrom+=1
			new_id=f"{record.id}[organism=Cryptococcus neoformans][strain={args.isolate_name}][location=chromosome][chromosome={chrom}]"
			new_seq=SeqRecord(record.seq,id=new_id,description='')
			SeqIO.write(new_seq,out_f,"fasta")

