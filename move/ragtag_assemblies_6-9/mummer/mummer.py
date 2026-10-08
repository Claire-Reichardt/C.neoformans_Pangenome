##Last updated:	11/10/25

##import system
import sys
import subprocess
import os
import argparse
import statistics
parser = argparse.ArgumentParser(description='Code to determine percent identity GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-n", metavar="NUCLEOTIDE_DIRECTORY", dest="nucdir", help="Directory containing all nucleotide assemblies to compare in .fasta files", required=True)
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
compdir = f"{userdir}/mummer/"

##Check files exist
NUC = args.nucdir
if os.path.isdir(NUC) and os.access(NUC, os.R_OK):
	genomedir = (os.path.abspath(args.nucdir))
else:
	sys.exit("ERROR: Assembly directory is missing or not readable!!")

##print files
genome_files = []
for fp in os.listdir(genomedir):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		genome_files.append(fp)
if len(genome_files)==0:
	sys.exit("".join(['ERROR: Unable to find any fasta files in: ', genomedir, '. ']))
else:
	print("".join(['Found ', str(len(genome_files)), ' asssembly files in: ', genomedir]))

##write shell files
total_id=[]
for file1 in os.listdir(genomedir):
	genome1=file1.split(".fa")[0]
	genome1_dest=f"{genomedir}/{file1}"
	mumdir=f"{compdir}{genome1}/"
	if not os.path.isdir(mumdir):
		shellfile=f"{mumdir}{genome1}_mshell.sh"
		os.makedirs(os.path.dirname(shellfile),exist_ok=True)
		shell = open(shellfile, 'w')
		shell.write ("""#!/bin/bash
#SBATCH --job-name=mummer
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=""")
		shell.write(mumdir)
		shell.write("""/%x_%j.out
""")
		shell.write("""#SBATCH --mail-user=""")
		if args.email is not None:
			shell.write(args.email)
			shell.write("\n")
			shell.write("""#SBATCH --mail-type=FAIL""")
		shell.write('\n')
		shell.write('\n')

		##write load modules
		shell.write("""module purge
ml MUMmer/4.0.0rc1-GCCcore-12.3.0
ml Python/3.11.3-GCCcore-12.3.0
ml Biopython/1.83-foss-2023a""")
		shell.write('\n')
		shell.write('\n')

		##write mummer functions
		shell.write("".join(['cd ', mumdir,'\n','\n']))
		for file2 in os.listdir(genomedir):
			if file1 != file2:
				genome2=file2.split(".fa")[0]
				genome2_dest=f"{genomedir}/{file2}"
				prefix=f"{genome1}_{genome2}"
				shell.write("".join(['nucmer -p ',prefix,' ',genome1_dest,
				' ',genome2_dest, '\n']))

				shell.write("".join(['dnadiff -p ',prefix,' -d ',
				prefix,'.delta','\n','\n']))

		shell.write("""find . -type f -not -name '*report' -print0 | xargs -0 rm --""")
		shell.write('\n')
		shell.write("".join(['python ',os.path.abspath(__file__),' -n ',genomedir,'\n']))
		##print to execute shell file
		shell.close()
		if args.submit is True:
			execute = ("".join(['sbatch ', shellfile, '\n']))
			os.system(execute)

	##get average percent id
	percent_id=[]
	for file in os.listdir(mumdir):
		ext=os.path.splitext(file)[-1].lower()
		report_file=f"{mumdir}{file}"
		if ext ==".report":
			ident_count=0
			with open(report_file) as fin:
				lines=fin.readlines()
			for line in lines:
				if "AvgIdentity" in line:
					ident_count+=1
					if ident_count==2:
						percid=float(line.split()[1])
						percent_id.append(percid)
	if len(percent_id)>0:
		average_percent_id=statistics.mean(percent_id)
		print("".join([genome1,"\t",str(average_percent_id)]))
		total_id.append(average_percent_id)

##Get percent id for all files
ave_total_id=statistics.mean(total_id)
ave_diverge=100-ave_total_id
print("".join(['Total percent id for assemblies is: ', str(ave_total_id), ' or percent divergence is: ',
str(ave_diverge)]))
