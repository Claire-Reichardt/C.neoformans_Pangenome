##Last updated:	10/1/24

##import system
import sys
import subprocess
import os
import argparse
from Bio import SeqIO
import pandas as pd
import scipy.stats as stats
parser = argparse.ArgumentParser(description='Code to get GO terms from annotations from annotations.py output on GACRC Slurm',
epilog='Written by Claire Reichardt')

##argument options
parser.add_argument('orthoannots', metavar='ORTHOLOG_ANNOTATIONS', help='ortho_annotations.tsv from annotations.py output')
parser.add_argument('classify', metavar='ORTHOLOG_CLASSIFICATION', help='classify.txt from classify.py output')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
core_enrich=f"{userdir}/core_goenrich.csv"
core_gocount=f"{userdir}/core_goorthocount.csv"
acc_enrich=f"{userdir}/accessory_goenrich.csv"
acc_gocount=f"{userdir}/accessory_goorthocount.csv"
sing_enrich=f"{userdir}/singleton_goenrich.csv"
sing_gocount=f"{userdir}/singleton_goorthocount.csv"

##Check files exist
ORTHOANNOTS = args.orthoannots
if os.path.isfile(ORTHOANNOTS) and os.access(ORTHOANNOTS, os.R_OK):
	orthoannots = (os.path.abspath(ORTHOANNOTS))
else:
	sys.exit("ERROR: ortho_annotations.tsv is missing or not readable!!")

CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.txt is missing or not readable!!")

#Get the classification for each ortholog
with open(classfile,'r') as fin:
	lines=fin.readlines()

class_dict={}
for line in lines:
	no_n=line.split("\n")[0]
	ortho=line[0:9]
	gene_class=no_n.split(",")[1]
	class_dict[ortho]=gene_class

#Go through each line in the .tsv
with open(orthoannots,'r') as fin2:
	lines2=fin2.readlines()

orthogo_dict={}
for line in lines2:
	ortho=line[0:9]
	annot=line.split("\t")[1]
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
gosing_dict={}
sing_count=0
for ortho,goterm in orthogo_dict.items():
	if ortho in class_dict:
		if class_dict[ortho]=="core":
			core_count+=1
			for term in goterm:
				if term not in gocore_dict:
					gocore_dict[term]=0
				gocore_dict[term]+=1
		if class_dict[ortho]=="accessory":
			acc_count+=1
			for term in goterm:
				if term not in goaccess_dict:
					goaccess_dict[term]=0
				goaccess_dict[term]+=1
		if class_dict[ortho]=="singleton":
			sing_count+=1
			for term in goterm:
				if term not in gosing_dict:
					gosing_dict[term]=0
				gosing_dict[term]+=1

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

	if goterm in gosing_dict:
		pos_sing_count=gosing_dict[goterm]
	else:
		pos_sing_count=0
	pos_non_core_count=pos_acc_count+pos_sing_count
	neg_non_core_count=acc_count+sing_count - pos_non_core_count

	#Do Fisher's test
	core_df=[[pos_core_count,pos_non_core_count],[neg_core_count,neg_non_core_count]]
	core_fisher_test=stats.fisher_exact(core_df, alternative='two-sided')
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

	if goterm in gosing_dict:
		pos_sing_count=gosing_dict[goterm]
	else:
		pos_sing_count=0
	pos_non_acc_count=pos_core_count+pos_sing_count
	neg_non_acc_count=core_count+sing_count - pos_non_acc_count

	#Do Fisher's test
	acc_df=[[pos_acc_count,pos_non_acc_count],[neg_acc_count,neg_non_acc_count]]
	acc_fisher_test=stats.fisher_exact(acc_df, alternative='two-sided')
	if pos_acc_count/neg_acc_count > pos_non_acc_count/neg_non_acc_count:
		acc_goenrich[goterm]=acc_fisher_test.pvalue

#Singleton
sing_goenrich={}
for goterm,count in gosing_dict.items():
	pos_sing_count=count
	neg_sing_count=sing_count-count

        #look for term in other categories
	if goterm in goaccess_dict:
		pos_acc_count=goaccess_dict[goterm]
	else:
		pos_acc_count=0

	if goterm in gocore_dict:
		pos_core_count=gocore_dict[goterm]
	else:
		pos_core_count=0
	pos_non_sing_count=pos_acc_count+pos_core_count
	neg_non_sing_count=acc_count+core_count - pos_non_sing_count

	#Do Fisher's test
	sing_df=[[pos_sing_count,pos_non_sing_count],[neg_sing_count,neg_non_sing_count]]
	sing_fisher_test=stats.fisher_exact(sing_df, alternative='two-sided')
	if pos_sing_count/neg_sing_count > pos_non_sing_count/neg_non_sing_count:
		sing_goenrich[goterm]=sing_fisher_test.pvalue


#print top ten values
core_top=dict(sorted(core_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])
acc_top=dict(sorted(acc_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])
sing_top=dict(sorted(sing_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])

#Print data to files
core_enrich_df=pd.DataFrame.from_dict(data=core_top, orient='index')
core_enrich_df.to_csv(core_enrich, sep=',')

acc_enrich_df=pd.DataFrame.from_dict(data=acc_top, orient='index')
acc_enrich_df.to_csv(acc_enrich, sep=',')

sing_enrich_df=pd.DataFrame.from_dict(data=sing_top, orient='index')
sing_enrich_df.to_csv(sing_enrich, sep=',')

core_gocount_df=pd.DataFrame.from_dict(data=gocore_dict, orient='index')
core_gocount_df.to_csv(core_gocount, sep=',')

acc_gocount_df=pd.DataFrame.from_dict(data=goaccess_dict, orient='index')
acc_gocount_df.to_csv(acc_gocount, sep=',')

sing_gocount_df=pd.DataFrame.from_dict(data=gosing_dict, orient='index')
sing_gocount_df.to_csv(sing_gocount, sep=',')
