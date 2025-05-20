#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import sys

def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		"extract sequences match with expression from plsdb", \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-i", action="store", dest="grep", 
					type=str, required=True, help="expression to grep between '\"'(REQUIRED)")
					
	parser.add_argument("-p", action="store", dest="PLSD_file", 
					type=str, required=True, help="PLSDB fasta file(REQUIRED)")
					
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
	
	PLSD_file = Arguments.PLSD_file
	fasta = open(PLSD_file,'r')
	lines = fasta.readlines()
	fasta.close()
	output = open(Arguments.output,'w')
	
	FLAG = False
	for line in lines :
		if line[0] == '>' :
			if Arguments.grep in line :
				FLAG = True
			else :
				FLAG = False
		if FLAG :
			output.write(line)
		
	output.close()


if __name__ == "__main__":
	main()	