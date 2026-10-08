##Last updated: 5/20/26

##import system
import sys
import subprocess
import os
import argparse
parser = argparse.ArgumentParser(description='Code that runs pilon and medaka polishing following canu for long read assemblies on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('--submit', help='Submit job directly to slurm', action='store_true')
parser.add_argument("-e", metavar="EMAIL", dest="email", help="Email for slurm output if job fails", default=None)
required = parser.add_argument_group("Required arguments", "You must provide all of the arguments listed below:")
required.add_argument('srr', metavar='SRAs', nargs='+', help='SRR names to download for individual strain')
required.add_argument("-n",dest='name', metavar='NAME', help='Name of isolate for finding previous canu data in current working directory',required=True)
required.add_argument("-l", metavar="LONG_FILE", dest="longfile", help="Long read fastq.gz file for medaka polishing", required=True)
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
assembly_dir=f"{userdir}/{args.name}/"
canu_dir=f"{assembly_dir}{args.name}_canu/"
canu_draft=f"{canu_dir}{args.name}.contigs.fasta"
medaka_dir=f"{assembly_dir}{args.name}_medaka/"
os.makedirs(os.path.dirname(medaka_dir), exist_ok=True)
medaka_draft=f"{medaka_dir}consensus.fasta"
pilon_dir=f"{assembly_dir}{args.name}_pilon/"
os.makedirs(os.path.dirname(pilon_dir), exist_ok=True)

##write shell file
if os.path.isfile(canu_draft):
	shellfile=f"{assembly_dir}pshell.sh"
	os.makedirs(os.path.dirname(shellfile),exist_ok=True)
	shell=open(shellfile,"w")
	shell.write("""#!/bin/bash
#SBATCH --job-name=pilon
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=3:00:00
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

	#Run medaka
	shell.write("""ml medaka/1.12.1-foss-2023a""")
	shell.write("".join(['\n','medaka_consensus -t 8 -i ',os.path.abspath(args.longfile),
		' -d ', canu_draft,' -o ', medaka_dir,'\n']))

	#Download fastqs of illumina SRAs
	shell.write("""module purge
ml SRA-Toolkit/3.0.3-gompi-2022a""")
	shell.write("".join(['\n', '\n']))

	shell.write("".join(['cd ',pilon_dir, '\n']))

	SRRlist = args.srr
	for srr in SRRlist:
		shell.write("".join(['fastq-dump --split-3 --gzip ',srr,'\n']))

	shell.write("".join(['cat *_1.fastq.gz > ', args.name, '_1.fq.gz','\n']))
	shell.write("".join(['cat *_2.fastq.gz > ', args.name, '_2.fq.gz','\n']))
	shell.write("".join(['rm *fastq.gz', '\n\n']))

	#loop through this 5x
	for i in range(1,6,1):

		shell.write("".join(['#Polish ',str(i),'\n']))
		#write load mapping modules
		shell.write("""module purge
ml SAMtools/1.16.1-GCC-11.3.0
ml minimap2/2.24-GCCcore-11.3.0""")
		shell.write("".join(['\n']))

		#write minimap
		id=f"{args.name}-{i}"
		sam=f"{id}_aln.sam"
		sam_sort=f"{id}_aln_sorted.sam"
		bam=f"{id}_aln.bam"
		if i==1:
			prev_draft=medaka_draft
		if i>1:
			prev_draft=f"{args.name}-{i-1}.fasta"
		new_draft=f"{args.name}-{i}"

		shell.write("".join(['minimap2 -ax sr ',prev_draft,' ',args.name,'_1.fq.gz ',
			args.name,'_2.fq.gz > ',sam,'\n']))
		shell.write("".join(['samtools sort ', sam, ' -o ',sam_sort, '\n']))
		shell.write("".join(['samtools view -b ', sam_sort,' > ', bam,'\n']))
		shell.write("".join(['samtools index ', bam,'\n']))

		#write pilon
		shell.write("""ml Pilon/1.24-Java-17
java -Xmx8G -jar $EBROOTPILON/pilon.jar""")
		shell.write("".join(['\n']))
		shell.write("".join(['pilon --genome ',prev_draft,' --frags ',bam,' --output ', new_draft, '\n']))
		shell.write("".join(['rm ',sam,' ',sam_sort,' ',bam,' ',bam,'.bai ','\n','\n','\n']))

		if i==5:
			shell.write("".join(['mv ',new_draft,".fasta ",f"{args.name}_polished.fasta",'\n']))

	##print to execute shell file
	shell.close()
	if args.submit is True:
		execute = ("".join(['sbatch ', shellfile, '\n']))
		os.system(execute)
