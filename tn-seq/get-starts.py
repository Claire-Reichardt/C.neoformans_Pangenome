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
parser = argparse.ArgumentParser(description='submodule of tn-seq.py',
epilog='Written by Claire Reichardt')

##arguments
parser.add_argument('dir')
args = parser.parse_args()

##file definitions
tn_dir = os.getcwd()

##Get final starts of each insertion from bam file
ref_bam=f"{tn_dir}/{args.dir}_ref.bam"
out_bed=f"{tn_dir}/{args.dir}.bed"
fwd_dict={}
rev_dict={}
if not os.path.isfile(out_bed):
	if os.path.isfile(ref_bam):
		samfile=pysam.AlignmentFile(ref_bam,"rb")
		reverse_reads=[read for read in samfile if read.is_reverse]
		for read in reverse_reads:
			loc=(read.reference_start)+1
			chrom=read.reference_name
			combo_name=f"{chrom}@{loc}+"
			read_name=read.query_name
			if combo_name in fwd_dict:
				fwd_dict[combo_name].append(read_name)
			else:
				fwd_dict[combo_name]=[read_name]
	if os.path.isfile(ref_bam):
		samfile=pysam.AlignmentFile(ref_bam,"rb")
		forward_reads=[read for read in samfile if not read.is_reverse]
		for read in forward_reads:
			loc = read.reference_end
			chrom=read.reference_name
			combo_name=f"{chrom}@{loc}-"
			read_name=read.query_name
			if combo_name in rev_dict:
				rev_dict[combo_name].append(read_name)
			else:
				rev_dict[combo_name]=[read_name]
	bed=open(out_bed,"w")
	for location, reads in fwd_dict.items():
		chrom=location.split("@")[0]
		loc=(location.split("@")[1]).split("+")[0]
		num_reads=str(len(reads))
		if chrom != 'None':
			bed.write("".join([chrom,'\t',loc,'\t',loc,'\t',(','.join(reads)),'\t',num_reads,'\t+\n']))
	for location, reads in rev_dict.items():
		chrom=location.split("@")[0]
		loc=(location.split("@")[1]).split("-")[0]
		num_reads=str(len(reads))
		if chrom != 'None':
			bed.write("".join([chrom,'\t',loc,'\t',loc,'\t',(','.join(reads)),'\t',num_reads,'\t-\n']))

	bed.close()

##Now get mapping statistics
data_list=[]
ds_bam=f"{tn_dir}/{args.dir}_sorted.bam"
ref_bam=f"{tn_dir}/{args.dir}_ref.bam"
all_reads=0
reads_with_ds=0
all_mapped_reads=0
reads_on_chr8=0
unique_mapping_locations=0
unique_sites_not_chr8=0
out_bed=f"{tn_dir}/{args.dir}.bed"
out_file=f"{tn_dir}/{args.dir}_stats.csv"
if os.path.isfile(out_bed) and os.stat(out_bed).st_size!=0:
	if os.path.isfile(ds_bam):
		ds_sam=pysam.AlignmentFile(ds_bam,"rb")
		for read in ds_sam:
			all_reads+=1
		ds_sam=pysam.AlignmentFile(ds_bam,"rb")
		mapped_reads=[read for read in ds_sam if not read.is_unmapped]
		for read in mapped_reads:
			reads_with_ds+=1
	if os.path.isfile(ref_bam):
		ref_sam=pysam.AlignmentFile(ref_bam,"rb")
		mapped_reads=[read for read in ref_sam if not read.is_unmapped]
		for read in mapped_reads:
			chrom=read.reference_name
			all_mapped_reads+=1
			if chrom=="Chr08":
				reads_on_chr8+=1
	with open(out_bed,"r") as fin:
		unique_sites=[]
		for line in fin:
			site=line.split("\t")[1]
			if site not in unique_sites:
				unique_sites.append(site)
				unique_mapping_locations+=1
				if "Chr08" not in line:
					unique_sites_not_chr8+=1
	#create final statistics
	per_reads_with_ds=reads_with_ds/all_reads
	per_ds_reads_mapped=all_mapped_reads/reads_with_ds
	per_sites_not_chr8=unique_sites_not_chr8/unique_mapping_locations
	per_mapping_chr8=reads_on_chr8/all_mapped_reads
	final_data=[args.dir,all_reads,reads_with_ds,per_reads_with_ds,all_mapped_reads,per_ds_reads_mapped,reads_on_chr8,per_mapping_chr8,unique_mapping_locations,unique_sites_not_chr8,per_sites_not_chr8]
	data_list.append(final_data)

#Create final dataframe
df=pd.DataFrame(final_data).T
df.columns=['Sample_name','Reads','Reads with ds','% Reads with ds','Reads with ds that map',
'% Reads with ds that map','Reads to Chr08','Percent mapping to Chr08','Unique mapping locations',
'Unique sites outside of Chr08','% of sites not on Chr08']

df.to_csv(out_file,index=False)
