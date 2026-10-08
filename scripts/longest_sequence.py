##Last updated: 1/17/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
parser = argparse.ArgumentParser(description='Code picks out the longest representative per orthogroup',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthodir', metavar='ORTHOGROUP_DIRECTORY', help='Orthogroup_Sequences folder from Orthofinder output')

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
longdir = f"{userdir}/longest_orthos.fa"
os.makedirs(os.path.dirname(longdir), exist_ok=True)

##Check files exist
ORTHODIR = args.orthodir
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(args.orthodir))
else:
	sys.exit("ERROR: Orthologous sequence directory is missing or not readable!!")

##Go through all the files
#All of the code below is taken directly from Mycopunk, see: https://github.com/MycoPunk/Afum_PopPan/blob/main/pipelines/pangenome/OrthoFinder/scripts/get_longest_protein_OrthoGroup.py

longest_seqs=[]
for file in os.listdir(orthodir):
	if file.endswith(".fa"):
		orthogroup = file.split(".")[0]
		filename = os.path.join(orthodir,file)
		longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
		for record in SeqIO.parse(filename, "fasta"):
			id=record.id
			genome=id.split("_")[0]
			if "CNAG" not in genome:
				if not longest or len(record) > len(longest):
					if genome not in longest.description:
						longest = record

		if longest:
			rename=SeqRecord(longest.seq, orthogroup, "", longest.id)
			longest_seqs.append(rename)

SeqIO.write(longest_seqs,longdir,'fasta')
