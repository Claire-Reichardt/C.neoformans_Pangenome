##Last updated:	8/8/25

##import system
import sys
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description='Script for converting blast output into a bed file on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
parser.add_argument('blastfile', metavar='BLAST_FILE', help='.txt File of outfmt 6 BLAST results')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()

##Check files exist
BLASTFILE = args.blastfile
if os.path.isfile(BLASTFILE) and os.access(BLASTFILE, os.R_OK):
	blastfile = (os.path.abspath(args.blastfile))
else:
	sys.exit("ERROR: BLAST input file is missing or not readable!!")

##Get all the information from each line of the blast output
filename=os.path.splitext(args.blastfile)[0]
bed_out=f"{userdir}/{filename}.bed"

with open(blastfile, 'r') as fin:
	lines=fin.readlines()
for line in lines:
	#blast names for outfmt 6
	qseqid=line.split("\t")[0]
	sseqid=line.split("\t")[1]
	pident=line.split("\t")[2]
	#length=line.split("\t")[3]
	#mismatch=line.split("\t")[4]
	#gapopen=line.split("\t")[5]
	qstart=line.split("\t")[6]
	qend=line.split("\t")[7]
	#sstart=line.split("\t")[8]
	#ssend=line.split("\t")[9]
	evalue=line.split("\t")[10]
	#bitscore=(line.split("\t")[11]).split("\n")[0]
	with open(bed_out, 'a') as f:
		f.write("".join([qseqid,'\t',qstart, '\t', qend, '\t', sseqid, '\t',evalue, '\n']))
fin.close()
f.close()
