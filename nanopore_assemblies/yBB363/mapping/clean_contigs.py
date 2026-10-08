##Last updated: 5/22/26

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
parser = argparse.ArgumentParser(description='Code to rename and reorient contigs to H99 reference sequence on GACRC Slurm',
epilog='Written by Claire Reichardt')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY_FILE", dest="asb_file", help="File of assembled draft contigs after polishing", required=True)
required.add_argument("-p",dest='paf_file', metavar='PAF_FILE', help='.paf file produced by miniprot2 after aligning draft genome to H99 refernce',required=True)
args = parser.parse_args()

#File definitons
userdir=os.getcwd()
isolate_name=os.path.splitext(args.asb_file)[0]
new_fasta=f"{userdir}/{isolate_name}_clean.fa"

##Define chromosome numbers
chroms={"CP003820.1":"Chr01",
	"CP003821.1":"Chr02",
	"CP003822.1":"Chr03",
	"CP003823.1":"Chr04",
	"CP003824.1":"Chr05",
	"CP003825.1":"Chr06",
	"CP003826.1":"Chr07",
	"CP003827.1":"Chr08",
	"CP003828.1":"Chr09",
	"CP003829.1":"Chr10",
	"CP003830.1":"Chr11",
	"CP003831.1":"Chr12",
	"CP003832.1":"Chr13",
	"CP003833.2":"Chr14"}

##Go thorugh draft sequences
asb_file=os.path.abspath(args.asb_file)
filtered_contigs=[]
if os.path.isfile(asb_file):
	for record in SeqIO.parse(asb_file,"fasta"):
		if len(record.seq)>200000:
			filtered_contigs.append(record.id)
	print("".join(['Found ', str(len(filtered_contigs)),' over 200,000bps in ', asb_file]))

##Go through paf file
paf_file=os.path.abspath(args.paf_file)
new_names=[]
new_seqs=[]
if os.path.isfile(paf_file):
	for contig_name in filtered_contigs:
		best_chrom=None
		best_length=0
		direction=None
		with open(paf_file,'r') as fin:
			for line in fin:
				if contig_name in line:
					contig_id=line.split("\t")[0]
					contig_length=int(line.split("\t")[1])
					strand=line.split("\t")[4]
					chrom_name=line.split("\t")[5]
					alignment_len=int(line.split("\t")[10])
					length=alignment_len/contig_length
					if length>best_length:
						act_chrom=chroms[chrom_name]
						best_length=length
						best_chrom=act_chrom
						direction=strand

		with open(new_fasta,'a') as out_f:
			for record in SeqIO.parse(asb_file,"fasta"):
				if best_chrom in new_names:
					best_chrom=f"{best_chrom}-2"
				if record.id==contig_name:
					if direction=="-":
						new_seq=record.reverse_complement(id=best_chrom,description='')
					if direction=="+":
						new_seq=SeqRecord(record.seq,id=best_chrom,description='')
					new_names.append(best_chrom)
					new_seqs.append(new_seq)
					SeqIO.write(new_seq,out_f,"fasta")

				#SeqIO.write(new_seqs,out_f,"fasta")
