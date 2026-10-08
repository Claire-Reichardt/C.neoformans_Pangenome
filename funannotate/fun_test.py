##Last updated:	1/5/25

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code runs batches of Funannotate on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')

parser.add_argument("-t", metavar="TRINITY", dest="trinity",
help="Reusing trinity.fasta file from previous run of 'Funannotate train' with the same RNA-seq reads.", default=None)

required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument("-a", metavar="ASSEMBLY", dest="assembly", help="Directory with assembly files", required=True)

required.add_argument("-f", metavar="FORWARD", dest="forward",
help="gzip file of forward RNA-seq reads", required=True)

required.add_argument("-r", metavar="REVERSE", dest="reverse",
help="gzip file of reverse RNA-seq reads", required=True)

required.add_argument("-e", metavar="EMAIL", dest="email", help="Email for remote searches", required=True)

required.add_argument("-bs", metavar="BUSCO_SPECIES", dest="buscospecies",
help="Species/phylum for BUSCO to use, see options with: funannotate outgroups -b", required=True)

required.add_argument("-as", metavar="AUGUSTUS_SPECIES", dest="augustusspecies",
help="Species/phylum for Augustus to use, see options with: funannotate species", required=True)

parser.add_argument("-p", metavar="PROTEIN_EVIDENCE", dest="proteinevidence",
help="Fasta file of known annotated proteins to predict proteins", default=None)

parser.add_argument("-te", metavar="TRANSCRIPT_EVIDENCE", dest="transcriptevidence",
help="Fasta file of mRNAs to align to predict proteins", default=None)

args = parser.parse_args()

##Check files exist
ASSPATH = args.assembly
if os.path.isdir(ASSPATH) and os.access(ASSPATH, os.R_OK):
	assembly = (os.path.abspath(args.assembly))
else:
	sys.exit("ERROR: assembly directory is missing or not readable!!")

FPATH = args.forward
if os.path.isfile(FPATH) and os.access(FPATH, os.R_OK):
	forward = (os.path.abspath(args.forward))
else:
	sys.exit("ERROR: Forward RNA-seq file is missing or not readable!!")

RPATH = args.reverse
if os.path.isfile(RPATH) and os.access(RPATH, os.R_OK):
	reverse = (os.path.abspath(args.reverse))
else:
	sys.exit("ERROR: Reverse RNA-seq file is missing or not readable!!")

##file definitions
userdir = os.getcwd()
jobsfile= f"{userdir}/jobs.txt"

##get files
count_files=0
strain_names=[]
for fp in os.listdir(assembly):
	ext = os.path.splitext(fp)[-1].lower()
	if ext == ".fa":
		count_files+=1
		strain = os.path.splitext(fp)[0]
		strain_names.append(strain)
if count_files==0:
	sys.exit("".join(['ERROR: Unable to find any assembly files in: ', assembly, '.']))
else:
	print("".join(['Found ', str(count_files), ' assembly files in: ', assembly]))

##Write shell files for every assembly
sorted_strains=sorted(strain_names)
count_shell = 0
job_numbers=[]

for strain in sorted_strains:
	count_shell+=1

	##file definitions
	header = strain[0:10].replace('_', '-')
	file=f"{strain}.fa"
	filename = os.path.join(assembly,file)
	sortfile = f"{strain}.sort.fasta"
	maskfile = f"{strain}.mask.fasta"
	fundir = f"{strain}_fun"
	straindir = f"{userdir}/{strain}"
	shellfile = f"{userdir}/{strain}/{strain}_fshell.sh"
	os.makedirs(os.path.dirname(shellfile), exist_ok=True)

	##create a file for each header
	shell = open(shellfile, "w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=3-00:00:00
#SBATCH --output=""")
	shell.write(straindir)
	shell.write("""/%x_%j.out""")
	shell.write("""
#SBATCH --mail-user=""")
	shell.write(args.email)
	shell.write("""
#SBATCH --mail-type=FAIL""")
	shell.write('\n')

	if count_shell >=25:
		if args.submit is True:
			shell.write("""#SBATCH --dependency=afterany:""")
			shell.write(str(job_numbers[0]))
			shell.write('\n')
			del job_numbers[0]

	shell.write('\n')
	shell.write('\n')

	##write load modules
	shell.write("""ml funannotate/1.8.15
source activate ${EBROOTFUNANNOTATE}""")
	shell.write('\n')
	shell.write('\n')

	##change directories
	shell.write("".join(['cd ', straindir, '\n', '\n']))

	##write funannotate scripts
	#sort
	shell.write("".join(['funannotate sort -i ', filename, ' -o ', sortfile, '\n']))

	#mask
	shell.write("".join(['funannotate mask -i ', sortfile, ' -o ', maskfile, '\n']))

	#train
	if args.trinity is not None:
		TRINITY=args.trinity
		if os.path.isfile(TRINITY) and os.access(TRINITY, os.R_OK):
			trinity = os.path.abspath(args.trinity)
			shell.write("".join(['ml DBD-mysql/4.050-GCC-11.3.0', '\n']))
			shell.write("".join(['funannotate train -i ', maskfile, ' -o ', fundir,
			' -l ', forward, ' -r ', reverse, ' --stranded RF --jaccard_clip --pasa_db mysql ',
			'--memory 40G --cpus 16 --trinity ', trinity, '\n']))
		else:
			sys.exit("ERROR: Trinity file is missing or not readable!!")
	else:
		shell.write("".join(['ml DBD-mysql/4.050-GCC-11.3.0', '\n']))
		shell.write("".join(['funannotate train -i ', maskfile, ' -o ', fundir,
		' -l ', forward, ' -r ', reverse, ' --stranded RF --jaccard_clip --pasa_db mysql ',
		'--memory 40G --cpus 16', '\n']))


	#check PASA ran correctly
	pasafile= f"{straindir}/{fundir}/training/funannotate_train.coordSorted.bam"
	shell.write("".join(['\n','if [[ -e ', pasafile, ' ]]; then', '\n', '\t',
	'echo "PASA ran correctly. Continuing."', '\n', 'else', '\n', '\t',
	'echo "CLAIRE ERROR: PASA has failed. Exiting now."', '\n', '\t', 'exit 1',
	'\n', 'fi', '\n', '\n']))

	#predict
	shell.write('\n')
	shell.write("""conda deactivate
module purge
ml funannotate/1.8.15
source activate ${EBROOTFUNANNOTATE}""")
	shell.write('\n')
	shell.write('\n')

	if args.proteinevidence is not None and args.transcriptevidence is not None:
		PROTEV=os.path.abspath(args.proteinevidence)
		TRANSEV=os.path.abspath(args.transcriptevidence)
		if os.path.isfile(PROTEV) and os.access(PROTEV, os.R_OK):
			if os.path.isfile(TRANSEV) and os.access(TRANSEV, os.R_OK):
				shell.write("".join(['funannotate predict -i ', maskfile, ' -o ', fundir,
				' -s ', args.augustusspecies,' --busco_db ', args.buscospecies,
				' --name ', header,' --strain ', strain, ' --protein_evidence ',PROTEV,
				' --augustus_species ', args.augustusspecies,
				' --transcript_evidence ', TRANSEV,' --cpus 24 --keep_no_stops', '\n']))

			else:
				sys.exit("ERROR: transcript evidence file is missing or not readable!!")
		else:
			sys.exit("ERROR: protein evidence file is missing or not readable!!")


	if args.proteinevidence is not None and args.transcriptevidence is None:
		PROTEV=os.path.abspath(args.proteinevidence)
		if os.path.isfile(PROTEV) and os.access(PROTEV, os.R_OK):
			shell.write("".join(['funannotate predict -i ', maskfile, ' -o ', fundir,
			' -s ', args.augustusspecies,' --busco_db ', args.buscospecies,
			' --name ', header,' --strain ', strain, ' --protein_evidence ',PROTEV,
			' --cpus 24 --keep_no_stops --augustus_species ', args.augustusspecies, '\n']))

		else:
			sys.exit("ERROR: protein evidence file is missing or not readable!!")

	if args.proteinevidence is None and args.transcriptevidence is not None:
		TRANSEV=os.path.abspath(args.transcriptevidence)
		if os.path.isfile(TRANSEV) and os.access(TRANSEV, os.R_OK):
			shell.write("".join(['funannotate predict -i ', maskfile, ' -o ', fundir,
			' -s ', args.augustusspecies,' --busco_db ', args.buscospecies,
			' --name ', header,' --strain ', strain, ' --transcript_evidence ',TRANSEV,
			' --cpus 24 --keep_no_stops --augustus_species ',args.augustusspecies,
			'\n']))
		else:
			sys.exit("ERROR: transcript evidence file is missing or not readable!!")

	if args.transcriptevidence is None and args.proteinevidence is None:
		shell.write("".join(['funannotate predict -i ', maskfile, ' -o ', fundir,
		' -s ', args.augustusspecies,' --busco_db ', args.buscospecies,
		' --name ', header,' --strain ', strain, ' --cpus 24 --min_protlen 35 --repeats2evm ',
		'--keep_no_stops --augustus_species ',args.augustusspecies,'\n']))

	#update
	#shell.write('\n')
	#shell.write("".join(['ml DBD-mysql/4.050-GCC-11.3.0', '\n']))
	#shell.write("".join(['funannotate update -i ', fundir, ' --cpus 24 --jaccard_clip ','\n']))

	#ipr
	shell.write('\n')
	shell.write("""conda deactivate
module purge
ml funannotate/1.8.15
source activate ${EBROOTFUNANNOTATE}""")
	shell.write('\n')
	shell.write('\n')

	shell.write("".join(['funannotate iprscan -i ', fundir,
	' -m local --iprscan_path /apps/eb/InterProScan/5.63-95.0-foss-2022a --cpus 24', '\n']))

	#run phobius independently
	home_folder = os.getenv('HOME')
	phobiuspath = f"{home_folder}/phobius.pl"
	if os.path.isfile(phobiuspath) and os.access(phobiuspath, os.R_OK):
		shell.write("".join([phobiuspath,' -short ', straindir,'/', fundir,'/predict_results/', args.augustusspecies,
		'_', strain, '.proteins.fa ', '> ', straindir, '/', fundir, '/predict_results/',
		'phobius.results.txt', '\n']))
	else:
		sys.exit(''.join(["ERROR: could not find phobius.pl in: ", phobiuspath, ". Please install Phobius directly and move phobius.pl executable to the correct location."]))

	#annotate
	shell.write("".join(['funannotate annotate -i ', fundir, ' --cpus 24 --busco_db ', args.buscospecies, ' --strain ', strain,
	' --phobius ', straindir, '/', fundir, '/predict_results/phobius.results.txt', '\n']))

	#remove subdirectories
	predict_dir=f"{straindir}/{fundir}/predict_misc/*/"
	train_dir=f"{straindir}/{fundir}/training/*/"
	update_dir=f"{straindir}/{fundir}/update_misc/*/"
	shell.write("".join(['rm -rf ',predict_dir, '\n']))
	shell.write("".join(['rm -rf ',train_dir, '\n']))
	#shell.write("".join(['rm -rf ',update_dir, '\n']))
	shell.close()

	#help with submitting:
	if args.submit is True:
		execute=("".join(['sbatch ', shellfile,'\n']))
		slurm_job_id=os.popen(execute).read()
		remove_sbatch=slurm_job_id.split("Submitted batch job ")[1]
		job_id=remove_sbatch.split("\n")[0]
		job_numbers.append(int(job_id))
		print("".join(['Submitted batch job ', job_id]))
