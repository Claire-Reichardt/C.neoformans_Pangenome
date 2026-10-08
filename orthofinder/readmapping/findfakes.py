##Last updated:	6/30/25

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import fnmatch
parser = argparse.ArgumentParser(description='Code to validate ortholog groups on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-s", metavar="SRRDIRECTORY", dest="srrdir", help="Directory containing both forward and reverse reads of all sequences formatted in .fq.gz", required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required.add_argument("-o",dest='orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)', required=True)
required.add_argument("-c",dest='cds', metavar='CD_SEQUENCES', help='a file containing all of the cds identified by funannotate',required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
missing=f"{userdir}/missing_orthos.csv"
longdir=f"{userdir}/longest_orthos.fa"
nucsfile = f"{userdir}/representative_nucletides.fa"
indexnucsfile=f"{userdir}/representative_nucletides.fa.amb"

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
SRR = args.srrdir
if os.path.isdir(SRR) and os.access(SRR, os.R_OK):
	srrdir = (os.path.abspath(args.srrdir))
else:
	sys.exit("ERROR: SRR directory is missing or not readable!!")

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

#Get total genomes used for analysis
countentries=0
for file in os.listdir(genomesdir):
	countentries +=1
total_genomes=str(countentries)
print("".join(['Found ', str(countentries), ' genomes in: ', genomesdir]))

#see if longestorthos has been run yet for these orthologs
if not os.path.isfile(longdir):
	print("Did not find longest missing ortholog file in this directory. Will run longest_sequence.py on missing sequences now")
	#get the protein sequences for the longest representative of these orthologs
	longest_seqs=[]
	for file in os.listdir(orthodir):
		if file.endswith(".fa"):
			orthogroup = file.split(".")[0]
			filename = os.path.join(orthodir,file)
			longest = SeqRecord(Seq("A"),id="placeholder",name="place",description="this is fake holder")
			for record in SeqIO.parse(filename, "fasta"):
				id=record.id
				genome=id.split("_")[0]
				protein=id.split("_")[1]
				if "CNAG" not in genome:
					if "A" not in protein or "B" not in protein:
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
					rename=SeqRecord(longest.seq, orthogroup, "", longest.id)
					longest_seqs.append(rename)
				if "CNAG" in id:
					id=longest.id
					new_id=id.split("-p")[0]
					rename=SeqRecord(longest.seq, orthogroup, "", new_id)
					longest_seqs.append(rename)

	SeqIO.write(longest_seqs,longdir,'fasta')
else:
	print("".join(['Found longest missing ortholog file: ', longdir]))

#Find representative nucleotides
if not os.path.isfile(nucsfile):
	#Get the longest sequence cd ID
	proteinIDs={}
	for record in SeqIO.parse(longdir, "fasta"):
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
else:
	print("".join(['Found representative nucleotide for missing ortholog file: ', nucsfile]))

#Write shell files for every strain to get read depth
count = 0
shellfiles=[]
for fp in os.listdir(srrdir):
	if fnmatch.fnmatch(fp, '*_1.fq.gz'):
		count +=1
		strain = fp.replace('_1.fq.gz','')
		forward = f"{srrdir}/{strain}_1.fq.gz"
		if os.path.isfile(forward) and os.access(forward, os.R_OK):
			forward = os.path.abspath(forward)
		reverse = f"{srrdir}/{strain}_2.fq.gz"
		if os.path.isfile(reverse) and os.access(reverse, os.R_OK):
			reverse= os.path.abspath(reverse)
		else:
			sys.exit(''.join(["ERROR:Could not find matching reverse reads for:",
			forward, ". Please make sure forward(_1) and reverse(_2) files are named the same."]))

		##file definitions
		shellfile = f"{userdir}/readmapping/{strain}/{strain}_map.sh"
		os.makedirs(os.path.dirname(shellfile), exist_ok=True)
		shellfiles.append(shellfile)

		straindir = f"{userdir}/readmapping/{strain}"
		forward_unzip=f"{straindir}/{strain}_1.fastq"
		reverse_unzip=f"{straindir}/{strain}_2.fastq"

		##create a file for each header
		shell = open(shellfile, "w")
		shell.write("""#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=""")
		shell.write(straindir)
		shell.write("""/%x_%j.out
""")
		if args.email is not None:
			shell.write("""#SBATCH --mail-user=""")
			shell.write(args.email)
			shell.write('\n')
			shell.write("""#SBATCH --mail-type=FAIL""")
		shell.write('\n')
		shell.write('\n')

		##write load modules
		shell.write("""module purge
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0""")
		shell.write('\n')
		shell.write('\n')

		##make index file if not yet done
		if not os.path.isfile(indexnucsfile):
			shell.write("".join(['bwa index ', nucsfile, '\n']))

		##unzip files
		shell.write("".join(['gunzip -ck ', forward, ' > ', forward_unzip, '\n']))
		shell.write("".join(['gunzip -ck ', reverse, ' > ', reverse_unzip, '\n']))

		##change directories
		shell.write("".join(['cd ', straindir, '\n']))
		shell.write('\n')

		##write commands
		shell.write("".join(['bwa mem ',nucsfile,' ',forward_unzip,' ',reverse_unzip,' > ',strain,'.sam','\n']))

		shell.write("".join(['samtools sort ',strain,'.sam -o ',strain,'_sorted.sam','\n' ]))

		shell.write("".join(['samtools view -b ',strain,'_sorted.sam > ',strain,'_sorted.bam','\n']))

		shell.write("".join(['samtools depth -a ',strain,'_sorted.sam -o ',strain, '_depth.sam','\n']))

		shell.write("".join(['samtools coverage ', strain, '_sorted.sam -o ', strain, '_coverage.tsv','\n']))

		shell.write("".join(['samtools view ', strain,
		'_sorted.bam | grep "XA" -> ', strain, '_multimap.sam','\n']))

		shell.close()

		##help with submission
		if args.submit is True:
			execute=("".join(['sbatch ', shellfile,'\n']))
			os.system(execute)
		else:
			print("".join(['--submit flag was not used. For submission of jobs to slurm please use: sbatch ',shellfile]))
