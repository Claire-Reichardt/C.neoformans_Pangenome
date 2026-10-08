##Last updated:	7/23/26

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import pandas as pd
import scipy.stats as stats
parser = argparse.ArgumentParser(description='Code to get GO terms from interproscan output on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('gff_file', metavar='GFF_FILE', help='gff file from interproscan output ')
parser.add_argument('classify', metavar='ORTHOLOG_CLASSIFICATION', help='classify.csv from new_classify.py output')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
core_enrich=f"{userdir}/core_goenrich.csv"
core_gocount=f"{userdir}/core_goorthocount.csv"
acc_enrich=f"{userdir}/accessory_goenrich.csv"
acc_gocount=f"{userdir}/accessory_goorthocount.csv"

##Check files exist
GFF_FILE = args.gff_file
if os.path.isfile(GFF_FILE) and os.access(GFF_FILE, os.R_OK):
	gfffile = (os.path.abspath(GFF_FILE))
else:
	sys.exit("ERROR: ortho_annotations.tsv is missing or not readable!!")

CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.csv is missing or not readable!!")

#Get the classification for each ortholog
class_dict={}
if os.path.isfile(classfile):
	with open(classfile,'r') as fin:
		lines=fin.readlines()
	for line in lines:
		no_n=line.split("\n")[0]
		ortho=line[0:9]
		gene_class=no_n.split(",")[1]
		class_dict[ortho]=gene_class

#Go through each line in the .tsv
orthogo_dict={}
if os.path.isfile(gfffile):
	with open(gfffile,'r') as fin:
		for line in fin:
			if "#" not in line and "\t" in line:
				ortho=line[0:9]
				annot=line.split("\t")[8]
				if "Ontology_term=" in annot:
					go_split=annot.split("Ontology_term=")[1]
					go_list=go_split.split(";")[0]
					go_terms=go_list.split(",")
					orthogo_dict[ortho]=list(go_terms)

#Put go terms in each classifer
gocore_dict={}
core_count=0
goaccess_dict={}
acc_count=0
for ortho,goterm in orthogo_dict.items():
	if ortho in class_dict:
		if class_dict[ortho]=="core":
			core_count+=1
			for term in goterm:
				term=term.replace('"','')
				if term not in gocore_dict:
					gocore_dict[term]=0
				gocore_dict[term]+=1
		if class_dict[ortho]=="accessory":
			acc_count+=1
			for term in goterm:
				term=term.replace('"','')
				if term not in goaccess_dict:
					goaccess_dict[term]=0
				goaccess_dict[term]+=1
print(core_count)

##Get enrichment numbers for each GOterm in each class
#core
core_goenrich={}
for goterm,count in gocore_dict.items():
	pos_core_count=count
	neg_core_count=core_count-count

	#look for term in other categories
	if goterm in goaccess_dict:
		pos_acc_count=goaccess_dict[goterm]
	else:
		pos_acc_count=0

	pos_non_core_count=pos_acc_count
	neg_non_core_count=acc_count - pos_non_core_count

	#Do Fisher's test
	core_df=[[pos_core_count,pos_non_core_count],[neg_core_count,neg_non_core_count]]
	core_fisher_test=stats.fisher_exact(core_df, alternative='two-sided')
	print("".join([str(pos_core_count), " ",str(pos_non_core_count), " ",str(neg_core_count)," ",
str(neg_non_core_count)," ",str(core_fisher_test.pvalue)]))
	if pos_core_count/neg_core_count > pos_non_core_count/neg_non_core_count:
		core_goenrich[goterm]=core_fisher_test.pvalue

#accessory
acc_goenrich={}
for goterm,count in goaccess_dict.items():
	pos_acc_count=count
	neg_acc_count=acc_count-count

	#look for term in other categories
	if goterm in gocore_dict:
		pos_core_count=gocore_dict[goterm]
	else:
		pos_core_count=0

	pos_non_acc_count=pos_core_count
	neg_non_acc_count=core_count - pos_non_acc_count

	#Do Fisher's test
	acc_df=[[pos_acc_count,pos_non_acc_count],[neg_acc_count,neg_non_acc_count]]
	acc_fisher_test=stats.fisher_exact(acc_df, alternative='two-sided')
	if pos_acc_count/neg_acc_count > pos_non_acc_count/neg_non_acc_count:
		acc_goenrich[goterm]=acc_fisher_test.pvalue


#print top ten values
core_top=dict(sorted(core_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])
acc_top=dict(sorted(acc_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])

#Print data to files
core_enrich_df=pd.DataFrame.from_dict(data=core_top, orient='index')
core_enrich_df.to_csv(core_enrich, sep=',')

acc_enrich_df=pd.DataFrame.from_dict(data=acc_top, orient='index')
acc_enrich_df.to_csv(acc_enrich, sep=',')

core_gocount_df=pd.DataFrame.from_dict(data=gocore_dict, orient='index')
core_gocount_df.to_csv(core_gocount, sep=',')

acc_gocount_df=pd.DataFrame.from_dict(data=goaccess_dict, orient='index')
acc_gocount_df.to_csv(acc_gocount, sep=',')
