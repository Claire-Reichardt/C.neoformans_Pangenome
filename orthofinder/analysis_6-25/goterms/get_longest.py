##Last updated: 7/22/26

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
parser = argparse.ArgumentParser(description='Code that runs miniprot to get longest protein representatives on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required.add_argument('-o',dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22/)', required=True)
args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
LOGFILE = f"{args.orthodir}/Log.txt"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
	logfile = (os.path.abspath(LOGFILE))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

##file definitions
userdir = os.getcwd()
longdir=f"{userdir}/longest_orthos.fa"

##print files
count = 0
for fp in os.listdir(orthodir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any Orthogroup files in: ', orthodir, '. Please check OrthoFinder ran completely in provided directory.']))
else:
	print("".join(['Found ', str(count), ' Orthogroup files in: ', orthodir]))

genomes=[]
headers=[]
with open(logfile, 'r') as fin:
	lines=fin.readlines()
for line in lines:
	if not "#" in line:
		if ":" in line:
			infile=line.split(": ")[1]
			ext=os.path.splitext(infile)[-1].lower()
			if ".fa" in ext or ".fasta" in ext:
				if "_" in infile:
					genomename=infile.split("_")[1]
				else:
					genomename=infile
				if genomename not in genomes:
					no_n=genomename.split("\n")[0]
					genomes.append(no_n)
					no_fa=os.path.splitext(no_n)[0]
					header=no_fa[0:10].replace('_','-')
					headers.append(header)
if len(genomes)>0:
	print("".join(['Found ', str(len(genomes)), ' genomes in : ', orthodir]))
else:
	sys.exit("ERROR: improperly formatted logfile in orthofinder directory")

##Get all orthologs
if not os.path.isfile(longdir):
	missing_orthos=[]
	for file in os.listdir(orthodir):
		genome_list=[]
		filename=("".join([orthodir,'/',file]))
		orthogroup=os.path.splitext(file)[0]
		missing_orthos.append(orthogroup)

	##Get representative sequence for each ortholog
	longest_seqs=[]
	for file in os.listdir(orthodir):
		orthogroup=os.path.splitext(file)[0]
		if orthogroup in missing_orthos:
			filename=("".join([orthodir,'/',file]))
			longest_non_CNAG=SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			longest_CNAG=SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			for record in SeqIO.parse(filename,"fasta"):
				id=record.id
				genome=id.split("_")[0]
				if "CNAG" in genome:
					if not longest_CNAG or len(record)> len(longest_CNAG):
						if genome not in longest_CNAG.description:
							longest_CNAG=record
				else:
					if not longest_non_CNAG or len(record)>len(longest_non_CNAG):
						if genome not in longest_non_CNAG.description:
							if "OG" not in record.description:
								if "_A" not in record.id and "_B" not in record.id:
									longest_non_CNAG=record
			if "CNAG" in longest_CNAG.description:
				rename_CNAG=SeqRecord(longest_CNAG.seq,orthogroup,"",longest_CNAG.id)
				longest_seqs.append(rename_CNAG)
			else:
				rename_non_CNAG=SeqRecord(longest_non_CNAG.seq,orthogroup,"",longest_non_CNAG.id)
				longest_seqs.append(rename_non_CNAG)

	SeqIO.write(longest_seqs,longdir,'fasta')
