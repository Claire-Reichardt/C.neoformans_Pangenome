##Last updated: 4/10/24

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that runs BUSCO for genome assemblies on GACRC Slurm',
epilog = 'Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument('--download', help ='Download busco lineage dataset if not in current directory', action='store_true')
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-s", metavar="STRAIN", dest="strain", help="Strain name for naming files", required=True)
required.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm progress", required=True)
required.add_argument("-a", metavar="ASSEMBLY", dest="assembly", help="Fasta genome assembly file", required=True)
required.add_argument("-l", metavar="LINEAGE", dest="lineage", help="Lineage dataset for busco to use, see options with: busco --list-datasets", required=True)

args = parser.parse_args()

##file definitions
userdir = os.getcwd()
strain = args.strain
lineagedir = f"{userdir}/busco_downloads/lineages/{args.lineage}"
straindir = f"{userdir}/{strain}"
shellfile = f"{userdir}/{strain}/{strain}_bshell.sh"
os.makedirs(os.path.dirname(shellfile), exist_ok=True)

##check assembly file is valid
ASSPATH = args.assembly
if os.path.isfile(ASSPATH) and os.access(ASSPATH, os.R_OK):
    sortfile= (os.path.abspath(args.assembly))
else:
    sys.exit("ERROR: Assembly file is missing or not readable!")
LPATH = lineagedir
if os.path.isdir(LPATH) and os.access(LPATH, os.R_OK):
    pass
elif args.download is True:
	pass
else:
    sys.exit("ERROR: BUSCO lineage folder does not exist in this directory. Please download with '--download'.")


##create a header for each file
shell = open(shellfile, 'w')
shell.write ("""#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --output=""")
shell.write(straindir)
shell.write("""/%x_%j.out
""")
shell.write("""#SBATCH --mail-user=""")
shell.write(args.email)
shell.write("""    # Where to send mail
#SBATCH --mail-type=FAIL"""
)
shell.write('\n')
shell.write('\n')

##write load modules
shell.write("""module purge
ml BUSCO/5.5.0-foss-2022a
usejni=t""")
shell.write('\n')
shell.write('\n')

##change directories
db = ['cd ', straindir, '\n']
printdb=''.join(db)
shell.write(printdb)
shell.write('\n')

##download BUSCO lineage dataset
if args.download is True:
        shell.write("".join(['busco --download ', args.lineage, '\n']))

##write BUSCO functions
run = ['busco -i ', sortfile, ' -m genome -l ', lineagedir, ' -o ', strain, '_busco -f', '\n']
printrun = ''.join(run)
shell.write(printrun)

##print to execute shell file
shell.close()
if args.submit is True:
	execute = ("".join(['sbatch ', shellfile, '\n']))
	os.system(execute)
