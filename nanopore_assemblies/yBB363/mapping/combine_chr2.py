##Last updated: 5/25/26

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
parser = argparse.ArgumentParser(description='Code to combine chrom2 with 50Ns after running clean_congits.fa on GACRC Slurm',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY_FILE", dest="asb_file", help="File of assembled draft contigs after clean_contigs.fa", required=True)
args = parser.parse_args()

#File definitons
userdir=os.getcwd()
isolate_name=os.path.splitext(args.asb_file)[0]
new_fasta=f"{userdir}/{isolate_name}_final.fa"

##Go through contigs
seqs=[]
asb_file=os.path.abspath(args.asb_file)
if os.path.isfile(asb_file):
	chr2_contigs=[]
	left_arm=None
	right_arm=None
	for record in SeqIO.parse(asb_file,"fasta"):
		if "Chr02" in record.id:
			chr2_contigs.append(record)
		else:
			seqs.append(record)
	for sequence in chr2_contigs:
		if len(sequence)> 500000:
			right_arm=sequence
		else:
			left_arm=sequence
	chr2_seq=Seq(left_arm.seq+"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"+right_arm.seq)
	chr2_record=SeqRecord(chr2_seq,id="Chr02",description='')
	seqs.append(chr2_record)
	sorted_new=sorted(seqs,key=lambda r:r.id)
	SeqIO.write(sorted_new,new_fasta,"fasta")
