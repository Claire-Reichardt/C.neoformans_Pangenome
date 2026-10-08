##Last updated:	12/10/25

##import system
import sys
import subprocess
import os
import argparse
import pysam
import pandas as pd
import fnmatch
import ast
parser = argparse.ArgumentParser(description='Code to map transposon insertions on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-t", metavar="TRANSPOSON_FILE", dest="tnfile", help="Fasta file containing transposon fragment", required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required.add_argument("-g",dest='genome', metavar='GENOME_FILE', help='Fasta file of reference sequence to map reads to', required=True)
required.add_argument("-s",dest='sequence', metavar='SEQUENCES_DIR', help='Directory of fastq sequence files of insertions, must be .fastq files',required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()

##Check files exist
SEQUENCE = args.sequence
if os.path.isdir(SEQUENCE) and os.access(SEQUENCE, os.R_OK):
	seqdir = (os.path.abspath(args.sequence))
else:
	sys.exit("ERROR: Sequence directory is missing or not readable!!")

GENOME = args.genome
if os.path.isfile(GENOME) and os.access(GENOME, os.R_OK):
	genomefile = (os.path.abspath(args.genome))
else:
	sys.exit("ERROR: genome file provided is missing or not readable!!")

TRANSPOSON = args.tnfile
if os.path.isfile(TRANSPOSON) and os.access(TRANSPOSON, os.R_OK):
        tnfile = (os.path.abspath(args.tnfile))
else:
        sys.exit("ERROR: transposon fragment file provided is missing or not readable!!")

##print files
count = 0
for fp in os.listdir(seqdir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fastq":
		count +=1
if count ==0:
	sys.exit("".join(['ERROR: Unable to find any .fastq files in: ', orthodir]))
else:
	print("".join(['Found ', str(count), ' sequence files in: ', seqdir]))

#Write shell files for every strain to get read depth
count = 0
shellfiles=[]

tn_dir=f"{userdir}/tnmapping/"

if not os.path.isdir(tn_dir):
	for fp in os.listdir(seqdir):
		if fnmatch.fnmatch(fp, '*.fastq'):
			count +=1
			id = fp.replace('.fastq','')

			shellfile = f"{tn_dir}{id}/{id}_map.sh"
			os.makedirs(os.path.dirname(shellfile), exist_ok=True)

			reads=f"{seqdir}/{fp}"
			iddir = f"{userdir}/tnmapping/{id}"

			##create a file for each header
			shell = open(shellfile, "w")
			shell.write("""#!/bin/bash
#SBATCH --job-name=tn_map
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=""")
			shell.write(iddir)
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
ml BWA/0.7.17-GCCcore-12.3.0
ml Pysam/0.22.0-GCC-12.3.0
ml Python/3.11.3-GCCcore-12.3.0
ml Biopython/1.83-foss-2023a
ml SAMtools/1.6-GCC-12.3.0
""")
			shell.write('\n')
			shell.write('\n')

			indextnfile=f"{tnfile}.amb"

			##make index file if not yet done
			if not os.path.isfile(indextnfile):
				if count==1:
					shell.write("".join(['bwa index ', tnfile, '\n']))
					print(fp)
			##change directories
			shell.write("".join(['cd ', iddir, '\n']))
			shell.write('\n')

			##write commands first tn mapping commands
			shell.write("".join(['bwa mem ',tnfile,' ',reads,' > ',id,'.sam','\n']))
			shell.write("".join(['samtools sort ',id,'.sam -o ',id,'_sorted.bam','\n' ]))
			shell.write("".join(['samtools view ',id,'_sorted.bam > ',id,'_sorted.sam','\n']))

			##execute softclip script
			shell.write("".join(['python ',userdir,'/soft-clip.py ',id,'\n']))

			##write second reference mapping
			directory=id
			full_dir=f"{tn_dir}{directory}"
			soft_clip=f"{full_dir}/{directory}_softclip.fastq"
			shell.write("".join(['bwa mem ', genomefile,' ',soft_clip,' > ',directory,'_ref.sam\n']))
			shell.write("".join(['samtools sort ',directory,'_ref.sam > ',directory,'_ref.bam\n']))
			shell.write("".join(['samtools index ',directory,'_ref.bam\n']))

			##execute get starts script
			shell.write("".join(['python ',userdir,'/get-starts.py ',id,'\n']))

			shell.close()

			##help with submission
			if args.submit is True:
				execute=("".join(['sbatch ', shellfile,'\n']))
				os.system(execute)
			else:
				print("".join(['--submit flag was not used. For submission of jobs to slurm',
				' please use: sbatch ',shellfile]))


