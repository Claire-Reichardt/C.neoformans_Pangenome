##Last updated: 4/23/24

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code to BLAST genome annotations to pangenome sequence GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm job failure", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('-p', metavar='PANGENOME', dest='pangenome',
help='Fasta file pangenome protein sequences for BLAST comparison', required=True)
required.add_argument('-r', metavar='REFERENCE',dest='reference',
help='Fasta file of reference genome protein sequences for BLAST comparison', required=True)
required.add_argument('-rh', metavar='REFERENCE_GENENAMES', dest='genename',
help='Gene naming system for reference (example=CNAG,Afu)', required=True)
args = parser.parse_args()

##Check files exist
PANGENOME = args.pangenome
if os.path.isfile(PANGENOME) and os.access(PANGENOME, os.R_OK):
	pangenome = (os.path.abspath(args.pangenome))
else:
	sys.exit("ERROR: Pangenome fasta sequence file is missing or not readable!!")
REFERENCE = args.reference
if os.path.isfile(REFERENCE) and os.access(REFERENCE, os.R_OK):
        reference = (os.path.abspath(args.reference))
else:
        sys.exit("ERROR: Reference genome fasta sequence file is missing or not readable!!")

##file definitions
userdir = os.getcwd()
blastdir = f"{userdir}/BLASTP/"
os.makedirs(os.path.dirname(blastdir), exist_ok=True)
shellfile = f"{blastdir}/blastshell.sh"
os.makedirs(os.path.dirname(shellfile), exist_ok=True)

##create a file for each header
shell = open(shellfile, "w")
shell.write("""#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=""")
shell.write(blastdir)
shell.write("""%x_%j.out
""")
if args.email is not None:
	shell.write("""#SBATCH --mail-user=""")
	shell.write(args.email)
	shell.write("""
#SBATCH --mail-type=FAIL""")
shell.write('\n')
shell.write('\n')

##write load modules
shell.write("""module purge
ml Python/3.10.4-GCCcore-11.3.0
ml BLAST+/2.13.0-gompi-2022a""")
shell.write("".join(['\n', '\n']))

##change directories
shell.write("".join(['cd ', blastdir, '\n', '\n']))

##make output file names
basename=os.path.basename(pangenome)
base=os.path.splitext(basename)[0]
out=f"{blastdir}{basename}"
output=f"{blastdir}{base}_blastP_results.txt"
outputfiltered =f"{blastdir}{base}_blastP_results_filtered.txt"
temp = f"{blastdir}temp.txt"
finalout = f"{blastdir}one_to_one_blast_hits.txt"

#Credits to Jason Stajich for the code below, see: https://github.com/MycoPunk/Afum_PopPan/blob/main/pipelines/pangenome/OrthoFinder/post-process/00_run_BLASTP.sh

##Run blast
#make blast db
shell.write("".join(['makeblastdb -in ', pangenome, ' -dbtype prot -out ',out, '\n']))

#Run blastp -note, we're setting the max_target_seqs to a redic high number here and filtering results afterwards to avoid known problems with the command line version of max_target_seqs ecxluding top hits. 
shell.write("".join(['blastp -query ', reference, ' -db ', out, ' -out ', output,
' -outfmt 6 -max_target_seqs 500 ',
'-soft_masking true -evalue 1e-05', '\n']))

#subset to get only the top hit - note these are already sorted in order of likelyhood (alignment length and eval)
#so we just want the lines where it's the first instance of the value in the first col (Af293 gene hits)
shell.write("".join(
["awk '!($1 in a){a[$1];print}' < ", output, " > ", outputfiltered, "\n"]))

#grab just the first two columns for analysis
shell.write("".join(
["cat ", outputfiltered, " | grep '^", args.genename,
"*' > ", temp, "\n"]))

shell.write("".join(
["cut -f 1,2 ", temp, " > ", finalout, "\n"]))

#remove temp file
shell.write("".join(['rm ', temp]))

##print to execute shell file
shell.close()
if args.submit is True:
	execute = ("".join(['sbatch ', shellfile, '\n']))
	os.system(execute)
