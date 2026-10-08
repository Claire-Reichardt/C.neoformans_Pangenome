##Last updated: 10/6/25

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that runs BUSCO for protein sequences and grabs their single copy orthologs on GACRC Slurm',
epilog = 'Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument('--download', help ='Download busco lineage dataset if not in current directory', action='store_true')

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm progress", default=None)
required.add_argument("-p", metavar="PROTEINS", dest="protein", help="Directory with fasta protein sequences", required=True)
required.add_argument("-l", metavar="LINEAGE", dest="lineage", help="Lineage dataset for BUSCO to use, see options with: busco --list-datasets", required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
lineagedir = f"{userdir}/busco_downloads/lineages/{args.lineage}"
SCO_dir=f"{userdir}/single_copy_orthologs/"
os.makedirs(os.path.dirname(SCO_dir), exist_ok=True)

##check assembly file is valid
PROTPATH = args.protein
if os.path.isdir(PROTPATH) and os.access(PROTPATH, os.R_OK):
    protpath= (os.path.abspath(args.protein))
else:
    sys.exit("ERROR: Protein file is missing or not readable!")
LPATH = lineagedir
if os.path.isdir(LPATH) and os.access(LPATH, os.R_OK):
    pass
elif args.download is True:
	pass
else:
    sys.exit("ERROR: BUSCO lineage folder does not exist in this directory. Please download with '--download'.")

#Create a shell for each protein
for file in os.listdir(protpath):
	filename=os.path.join(protpath,file)
	if "cryptococcus" in file:
		id=file.split("cryptococcus_")[1]
		strain=id.split(".proteins.fa")[0]
	else:
		if ".fa" in file:
			strain=file.split(".fa")[0]
	straindir = f"{userdir}/{strain}"
	shellfile = f"{userdir}/{strain}/{strain}_bshell.sh"
	os.makedirs(os.path.dirname(shellfile), exist_ok=True)

	shell = open(shellfile, 'w')
	shell.write ("""#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=""")
	shell.write(straindir)
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
ml BUSCO/5.8.3-foss-2023a
usejni=t""")
	shell.write('\n')
	shell.write('\n')

	##change directories
	shell.write("".join(['cd ', straindir, '\n', '\n']))

	##download BUSCO lineage dataset
	if args.download is True:
 	       shell.write("".join(['busco --download ', args.lineage, '\n']))

	##write BUSCO functions
	shell.write("".join(['busco -i ', filename, ' -m prot -l ', lineagedir, ' -o ', strain, '_busco -f -c 8',
	'\n']))

	##combine single copy orthologs at the end of BUSCO run
	shell.write("".join(['cat ',strain,
	'_busco/run_basidiomycota_odb12/busco_sequences/single_copy_busco_sequences/* > ',
	strain,'_SCO.fa', '\n']))

	##move files
	shell.write("".join(['mv ',strain,'_SCO.fa ', SCO_dir, '\n']))

	##print to execute shell file
	shell.close()
	SCO_file=f"{SCO_dir}{strain}_SCO.fa"
	if not os.path.isfile(SCO_file):
		if args.submit is True:
			execute = ("".join(['sbatch ', shellfile, '\n']))
			os.system(execute)
