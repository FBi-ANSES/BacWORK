#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import os, sys
import json


def get_parser() :

	parser = argparse.ArgumentParser(description= \
		"convert tsv file to config.file for ARTwork2 snakemake", \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-i", action="store", dest="tab", 
					type=str, required=True, help="tab file(REQUIRED)")
					
	parser.add_argument("-dir", action="store", dest="reads_dir", 
					type=str, required=True, help="reads directory(REQUIRED)")
					
	parser.add_argument("-o", action="store", dest="job_dir", 
					type=str, required=True, help="working directory(REQUIRED)")
	
	return parser


def main():


	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
		
	Arguments=parser.parse_args()
	
	json_dico = {}
	json_dico["WORKDIR"] = os.path.abspath(Arguments.job_dir)
	json_dico["FASTQ_DIR"] = Arguments.reads_dir
	json_dico["samples"] = []
	
	try:
		f = open(Arguments.tab,'r')
		lines = f.readlines()
		f.close()
		for line in lines :
			line = line.rstrip().split()
			if ("SampleID" in line) or ("Sequencing_center" in line) or ("Genus" in line):
				continue
			else :
				dico_sample = {}
				dico_sample["SampleID"] = line[0]
				dico_sample["Project"] = line[1]
				dico_sample["Supplier"] = line[2]
				dico_sample["Sequencing_center"] = line[3]
				dico_sample["Phylogeny"] = {}
				dico_sample["Phylogeny"]["Genus"] = line[4]
				dico_sample["Phylogeny"]["Species"] = line[5]
				json_dico["samples"].append(dico_sample)
	except:
		exit()
		
	with open(f"config.json", "w") as outfile:  
		json.dump(json_dico, outfile, indent = 4)
				
				
				
	
	
if __name__ == "__main__":
	main()		