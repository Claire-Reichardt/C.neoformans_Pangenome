##Last updated:	8/19/24

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
parser.add_argument('orthodir', metavar='ORTHFINDER_DIRECTORY', help='Orthofinder output directory (ex.Orthofinder_diamond/Results_Apr22_1/)')
parser.add_argument('orthoannots', metavar='ORTHOLOG_ANNOTATIONS', help='ortho_annotations.tsv from annotations.py output')
parser.add_argument('classify', metavar='ORTHOLOG_CLASSIFICATION', help='classify.txt from classify.py output')
parser.add_argument('envkey', metavar='ENVIRONMENT_KEY', help='.CSV file containing header used for proteins of each genome and the clinical/environmental isolation associated with it')
args = parser.parse_args()

##file definitions
userdir = os.getcwd()
clin_enrich=f"{userdir}/clin_goenrich.csv"
clin_gocount=f"{userdir}/clin_goorthocount.csv"
env_enrich=f"{userdir}/env_goenrich.csv"
env_gocount=f"{userdir}/env_goorthocount.csv"

##Check files exist
ORTHODIR = f"{args.orthodir}/Orthogroup_Sequences"
GENOMESDIR =f"{args.orthodir}/all_orthos"
ENVKEY=args.envkey
if os.path.isdir(ORTHODIR) and os.access(ORTHODIR, os.R_OK):
	orthodir = (os.path.abspath(ORTHODIR))
else:
	sys.exit("ERROR: Orthogroup sequence directory is missing or not readable!!")

if os.path.isdir(GENOMESDIR) and os.access(GENOMESDIR, os.R_OK):
	genomesdir = (os.path.abspath(GENOMESDIR))
else:
	sys.exit("ERROR: getallorthos.py was not run for this Orthofinder run. Please execute getallorthos.py to ensure singleton genes are included in gene counts!!")
if os.path.isfile(ENVKEY) and os.access(ENVKEY, os.R_OK):
	envkey = (os.path.abspath(ENVKEY))
else:
	sys.exit("ERROR: matkey.csv provided is missing or not readable!!")
CLASSIFY = args.classify
if os.path.isfile(CLASSIFY) and os.access(CLASSIFY, os.R_OK):
	classfile = (os.path.abspath(CLASSIFY))
else:
	sys.exit("ERROR: classify.txt is missing or not readable!!")

ORTHOANNOTS = args.orthoannots
if os.path.isfile(ORTHOANNOTS) and os.access(ORTHOANNOTS, os.R_OK):
	orthoannots = (os.path.abspath(ORTHOANNOTS))
else:
	sys.exit("ERROR: ortho_annotations.tsv is missing or not readable!!")

#Get the classification for each ortholog
with open(classfile,'r') as fin:
	lines=fin.readlines()

class_dict={}
for line in lines:
	no_n=line.split("\n")[0]
	ortho=line[0:9]
	gene_class=no_n.split(",")[1]
	class_dict[ortho]=gene_class

#Read input csv to get key
env_dict={}
with open(envkey, 'r') as f:
	for line in f:
		items=line.split(',')
		key, values = items[0], items[1]
		values_2= values.split('\n')[0]
		env_dict[key]=values_2

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

#Get the orthogs for each genome
genome_ortholog_dict={}
for file in os.listdir(genomesdir):
	genomename=file.split("_")[0]
	name=os.path.splitext(genomename)[0]
	filename=f"{genomesdir}/{file}"

	unique_orthos=[]
	if "CNAG" not in filename:
		with open(filename, 'r') as fin:
			lines=fin.readlines()
		for line in lines:
			ortho=line[0:9]
			if ortho not in unique_orthos:
				unique_orthos.append(ortho)
	genome_ortholog_dict[name]=unique_orthos

#Put go terms in each classifer
goclin_dict={}
unique_orthos_clin=[]
clin_count=0
goenv_dict={}
unique_orthos_env=[]
env_count=0

for genome,ortho in genome_ortholog_dict.items():
	if genome in env_dict:
		if env_dict[genome]=="Clinical":
			for ortholog in ortho:
				if ortholog in orthogo_dict.keys():
					if ortholog not in unique_orthos_clin:
						unique_orthos_clin.append(ortholog)
						clin_count+=1
						goterms=list(orthogo_dict[ortholog])
						for term in goterms:
							if term not in goclin_dict:
								goclin_dict[term]=0
							goclin_dict[term]+=1
		if env_dict[genome]=="Environmental":
			for ortholog in ortho:
				if ortholog in orthogo_dict.keys():
					if ortholog not in unique_orthos_env:
						unique_orthos_env.append(ortholog)
						env_count+=1
						goterms=list(orthogo_dict[ortholog])
						for term in goterms:
							if term not in goenv_dict:
								goenv_dict[term]=0
							goenv_dict[term]+=1

##Get enrichment numbers for each GOterm in each class
#clinical
clin_goenrich={}
for goterm,count in goclin_dict.items():
	pos_clin_count=count
	neg_clin_count=clin_count-count

	#look for term in other categories
	if goterm in goenv_dict:
		pos_env_count=goenv_dict[goterm]
	else:
		pos_env_count=0

	pos_non_clin_count=pos_env_count
	neg_non_clin_count=env_count - pos_non_clin_count

	#Do Fisher's test
	clin_df=[[pos_clin_count,pos_non_clin_count],[neg_clin_count,neg_non_clin_count]]
	clin_fisher_test=stats.fisher_exact(clin_df, alternative='greater')
	if pos_clin_count/neg_clin_count > pos_non_clin_count/neg_non_clin_count:
		clin_goenrich[goterm]=clin_fisher_test.pvalue
		print(clin_df)
		print(clin_fisher_test.pvalue)

#environmental
env_goenrich={}
for goterm,count in goenv_dict.items():
	pos_env_count=count
	neg_env_count=env_count-count

	#look for term in other categories
	if goterm in goclin_dict:
		pos_clin_count=goclin_dict[goterm]
	else:
		pos_clin_count=0

	pos_non_env_count=pos_clin_count
	neg_non_env_count=clin_count - pos_non_env_count

	#Do Fisher's test
	env_df=[[pos_env_count,pos_non_env_count],[neg_env_count,neg_non_env_count]]
	env_fisher_test=stats.fisher_exact(env_df, alternative='greater')
	if pos_env_count/neg_env_count > pos_non_env_count/neg_non_env_count:
		env_goenrich[goterm]=env_fisher_test.pvalue


#print top ten values
clin_top=dict(sorted(clin_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])
env_top=dict(sorted(env_goenrich.items(), key = lambda x: x[1], reverse = False)[:10])

#Print data to files
clin_enrich_df=pd.DataFrame.from_dict(data=clin_top, orient='index')
clin_enrich_df.to_csv(clin_enrich, sep=',')

env_enrich_df=pd.DataFrame.from_dict(data=env_top, orient='index')
env_enrich_df.to_csv(env_enrich, sep=',')

clin_gocount_df=pd.DataFrame.from_dict(data=goclin_dict, orient='index')
clin_gocount_df.to_csv(clin_gocount, sep=',')

env_gocount_df=pd.DataFrame.from_dict(data=goenv_dict, orient='index')
env_gocount_df.to_csv(env_gocount, sep=',')
