#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import sys

def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		"remove plasmid from chromosome fasta", \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-c", action="store", dest="chr", 
					type=str, required=True, help="chromosome fasta(REQUIRED)")
					
	parser.add_argument("-p", action="store", dest="pls", 
					type=str, required=True, help="plasmid fasta file(REQUIRED)")
					
	parser.add_argument("-o", action="store", dest="output", 
					type=str, default="output", help="output file name (default:output)")

	return parser


def main():

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
		
	Arguments=parser.parse_args()
	
	file = Arguments.pls
	fasta = open(file,'r')
	lines = fasta.readlines()
	fasta.close()
	output = open(Arguments.output,'w')
	header = []
	
	for line in lines :
		if line[0] == '>' :
			header.append(line)
	
	file = Arguments.chr
	fasta = open(file,'r')
	lines = fasta.readlines()
	fasta.close()
	
	FLAG = True
	for line in lines :
		if line[0] == '>' :
			if line in header :
				FLAG = False
			else :
				FLAG = True
				header.append(line)
		if FLAG :
			output.write(line)
		
	output.close()


if __name__ == "__main__":
	main()	