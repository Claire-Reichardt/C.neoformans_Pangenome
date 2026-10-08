##Last updated:	7/1/24

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code average read depth for a contig for output from samtools depth on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('samdepth', metavar='SAMTOOLS_DEPTH_FILE', help='File from ouput of samtools depth')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
outfile = f"{userdir}/samtools_average_depth.tsv"

##Check files exist
SAMDEPTH = args.samdepth
if os.path.isfile(SAMDEPTH) and os.access(SAMDEPTH, os.R_OK):
	samfile = (os.path.abspath(SAMDEPTH))
else:
	sys.exit("ERROR: Samtools depth file is missing or not readable!!")

##Get the average read depth per ortholog
ortho_counter={}
read_counter={}

with open(samfile,'r') as fin:
	lines=fin.readlines()

for line in lines:
	ortho=line[0:9]
	depth=int(line.split("\t")[2])
	if ortho not in ortho_counter:
		ortho_counter[ortho]=depth
		read_counter[ortho]=0
	read_counter[ortho]+=1

	if ortho in ortho_counter:
		ortho_counter[ortho]+=(depth)

for ortho in ortho_counter:
	sum=ortho_counter[ortho]
	number=read_counter[ortho]
	average=sum/number

	with open(outfile, 'a') as f:
		f.write("".join([ortho, '\t', str(average), '\n']))

