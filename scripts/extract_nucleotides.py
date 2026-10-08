##Last updated: 7/1/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
parser = argparse.ArgumentParser(description='Code to grab nucleotide sequences from longest_sequences.py output and cds',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('longest_proteins', metavar='LONGEST_PROTEINS', help='longest_orthos.fa output from longest_sequences.py')
parser.add_argument('cds', metavar='CD_SEQUENCES', help='a file containing all of the cds identified by funannotate')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
nucsfile = f"{userdir}/representative_nucletides.fa"
os.makedirs(os.path.dirname(nucsfile), exist_ok=True)

##Check files exist
PROTFILE = args.longest_proteins
if os.path.isfile(PROTFILE) and os.access(PROTFILE, os.R_OK):
	protfile = (os.path.abspath(args.longest_proteins))
else:
	sys.exit("ERROR: longest_proteins file provided is missing or not readable!!")

CDS = args.cds
if os.path.isfile(CDS) and os.access(CDS, os.R_OK):
	cdsfile = (os.path.abspath(args.cds))
else:
	sys.exit("ERROR: cds file provided is missing or not readable!!")

#Get the longest sequence cd ID
proteinIDs={}
for record in SeqIO.parse(protfile, "fasta"):
	full=record.description
	proteinID=full.split(" ")[1]
	ortho=full.split(" ")[0]
	proteinIDs[proteinID]=ortho

nuc_seqs=[]
for record in SeqIO.parse(cdsfile, "fasta"):
	if record.id in proteinIDs:
		ortho=proteinIDs[record.id]
		rename=SeqRecord(record.seq, ortho, "", record.id)
		nuc_seqs.append(rename)

SeqIO.write(nuc_seqs,nucsfile,'fasta')
