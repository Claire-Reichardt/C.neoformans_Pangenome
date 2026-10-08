##Last updated: 5/20/26

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that runs canu for long read assemblies on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-l", metavar="LONG_FILE", dest="longfile", help="Long read fastq file for canu assembly", required=True)
required.add_argument("-n",dest='name', metavar='NAME', help='Name of isolate for prefixes',required=True)
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
assembly_dir=f"{userdir}/{args.name}/"

##write shell file
if os.path.isfile(args.longfile):
	long_dir=os.path.abspath(args.longfile)
	shellfile=f"{assembly_dir}cshell.sh"
	os.makedirs(os.path.dirname(shellfile),exist_ok=True)
	shell=open(shellfile,"w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=""")
	shell.write(f"{assembly_dir}")
	shell.write("""%x_%j.out
""")
	if args.email is not None:
		shell.write("""#SBATCH --mail-user=""")
		shell.write(args.email)
		shell.write("""
#SBATCH --mail-type=FAIL""")
	shell.write('\n')
	shell.write('\n')

	#write load modules
	shell.write("""module purge
cd $SLURM_SUBMIT_DIR
ml canu/2.3-GCCcore-13.3.0-Java-11""")
	shell.write("".join(['\n', '\n']))

	#write canu executions
	shell.write(''.join(['canu gridOptions=',
		'" --partition=batch --ntasks=1 --cpus-per-task=4 ',
		'--time=168:00:00 --mem-per-cpu=8G" -p ', args.name,
		' -d ',f"{assembly_dir}{args.name}_canu", ' contigFilter="2 200000 1.0 0.5 5"',
		' correctedErrorRate=0.075 genomeSize=19m -nanopore ',
		long_dir,'\n']))

	shell.close()
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)
