##Last updated:	12/10/25

##import system
import sys
import subprocess
import os
import pysam
import argparse
import pandas as pd
import fnmatch
import ast
parser = argparse.ArgumentParser(description='Submodule of tn-map.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
args = parser.parse_args()

##file definitions
dir = os.getcwd()

##extract first soft clip to a fastq file
soft_clip=f"{dir}/{args.dir}_softclip.fastq"
soft_file=open(soft_clip,"w")
for file in os.listdir(dir):
	if "_sorted.sam" in file:
		sam_sort=f"{dir}/{file}"
		with open(sam_sort,'r') as fin:
			for line in fin:
				name=line.split("\t")[0]
				cigar=line.split("\t")[5]
				seq=line.split("\t")[9]
				qual=line.split("\t")[10]
				if cigar != "*":
					if "S" in cigar:
						s_string=int(cigar.split("S")[0])
						soft=seq[0:s_string]
						soft_qual=qual[0:s_string]
						soft_file.write("".join(["@",name,'\n',soft,
						'\n','+','\n',soft_qual,'\n']))


