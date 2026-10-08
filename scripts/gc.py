##Last updated: 9/20/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from statistics import mean
from pathlib import Path
import pandas as pd
parser = argparse.ArgumentParser(description='Code to calculate GC content of orthologs in comparison to nanopore long read sequences on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-o",dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)', required=True)
required.add_argument("-c",dest='cds', metavar='CD_SEQUENCES', help='a file containing all of the cds identified by funannotate',required=True)
required.add_argument("-n",dest='nano', metavar='NANOPORE_FILE', help='a nucleotide file of long read nanopore sequences to compare GC content of orthologs to',required=True)

args = parser.parse_args()

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
GENOMESDIR =f"{args.orthodir}/all_orthos"
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")
CDS = args.cds
if os.path.isfile(CDS) and os.access(CDS, os.R_OK):
	cdsfile = (os.path.abspath(args.cds))
else:
	sys.exit("ERROR: cds file provided is missing or not readable!!")
NANO = args.nano
if os.path.isfile(NANO) and os.access(NANO, os.R_OK):
        nanofile = (os.path.abspath(args.nano))
else:
        sys.exit("ERROR: nanopore file provided is missing or not readable!!")

##file definitions
userdir = os.getcwd()
genome=Path(nanofile).stem
nucsfile=f"{userdir}/{genome}_repnucs.fa"
gc_file=f"{userdir}/{genome}_gc.csv"

##Calculate GC of nanopore genome
gc_list=[]
for record in SeqIO.parse(nanofile,"fasta"):
	gc=GC(record.seq)
	#remove mitochondrial genomes
	if gc > 40:
		gc_list.append(gc)
ave_gc=mean(gc_list)
print("".join(['Found ', str(len(gc_list)), ' contigs in ', nanofile, ' with an average GC of ', str(ave_gc)]))

#Find orthologs for the nanopore sequence
orthos_tofind=[]
proteinIDs={}
for file in os.listdir(genomesdir):
	if genome in file:
		filename = os.path.join(genomesdir,file)
		with open(filename,'r') as fin:
			lines=fin.readlines()
		for line in lines:
			no_n=line.split("\n")[0]
			ortho=line[0:9]
			if "\t" in line:
				if "CNAG" not in line:
					protein=no_n.split("\t")[1]
					proteinIDs[protein]=ortho
				if "CNAG" in line:
					protein=no_n.split("\t")[1]
					new_id=protein.split("-p")[0]
					proteinIDs[new_id]=ortho
			else:
				orthos_tofind.append(ortho)

#Find longest orthologs for orthos added
for ortho in orthos_tofind:
	filename=os.path.join(orthodir,f"{ortho}.fa")
	if os.path.isfile(filename):
		longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
		for record in SeqIO.parse(filename, "fasta"):
			id=record.id
			genome=id.split("_")[0]
			if "CNAG" not in genome:
				if not longest or len(record) > len(longest):
					if genome not in longest.description:
						longest = record
			if "CNAG" in genome:
				cnag_genome=genome.split("-p")[0]
				if not longest or len(record) > len(longest):
					if cnag_genome not in longest.description:
						longest = record

		if longest:
			if "CNAG" not in id:
				rename=SeqRecord(longest.seq, ortho, "", record.id)
				proteinIDs[record.id]=ortho
			if "CNAG" in id:
				id=record.id
				new_id=id.split("-p")[0]
				rename=SeqRecord(longest.seq, ortho, "", new_id)
				proteinIDs[new_id]=ortho

#Find nucleotide sequences for each protein
nuc_seqs=[]
for record in SeqIO.parse(cdsfile, "fasta"):
	if record.id in proteinIDs:
		ortho=proteinIDs[record.id]
		rename=SeqRecord(record.seq, ortho, "", record.id)
		nuc_seqs.append(rename)

SeqIO.write(nuc_seqs,nucsfile,'fasta')

#Compare GC of nucs to nanopore
gc_dict={}
bad_gcs=[]
for record in SeqIO.parse(nucsfile, "fasta"):
	gc=GC(record.seq)
	upper_gc=ave_gc+5
	lower_gc=ave_gc-5
	gc_dict[record.id]=gc-ave_gc
	if gc>upper_gc or gc<lower_gc:
		bad_gcs.append(record)

#Print GCs to file
gc_df=pd.DataFrame.from_dict(data=gc_dict, orient='index')
gc_df.to_csv(gc_file, sep=',')
