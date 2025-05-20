#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import json, os, re, sys


def get_parser() :

	parser = argparse.ArgumentParser(description= \
		"parse confindr and kraken2 output to detect contamination", \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-c", action="store", dest="confindr", 
					type=str, required=True, help="confindr report file(REQUIRED)")
					
	parser.add_argument("-k", action="store", dest="kraken2", 
					type=str, required=True, help="kraken2 output(REQUIRED)")
					
	parser.add_argument("-g", action="store", dest="genus", 
					type=str, required=True, help="genus(REQUIRED)")

	parser.add_argument("-s", action="store", dest="species", 
					type=str, required=True, help="species(REQUIRED)")					
					
	parser.add_argument("-t", action="store", dest="confindr_SNVs_threshold", 
					type=int, required=True, help="confindr SNVs threshold (REQUIRED)")
					
	parser.add_argument("-o", action="store", dest="output", 
					type=str, required=True, help="output file (REQUIRED)")				
	
	return parser



#Load confindr results
#Contamination by confindr is reported if both R1 and R2 are found contaminated
def Confindr_contamination(f,confindr_SNVs_threshold):
	with open(f) as f:
		header = f.readline().rstrip()
		line_R1 = f.readline().rstrip().split(',')
		confindr = 'False'

		if int(line_R1[2]) > confindr_SNVs_threshold:
			confindr = 'True'

	return(confindr)
   
#Load kraken results
#BAD quality kraken if percent secondary genus or species > 2%
def kraken_extract(f, threshold_pct):
	"""
	Report all genus and species with a percent superior to 2%
	"""

	with open(f) as kraken:
		genus_conta, species_conta = [], []
		for line in kraken:
			line = line.rstrip().split('\t')
			if line[3] == 'G':
				pct = round(float(line[0]),2)
				if pct > threshold_pct:
					species = line[-1]
					species = re.search('([A-Z].*)', species).group(1)
					to_report = species+':'+str(pct)
					genus_conta.append(to_report)

			elif line[3] == 'S':
				pct = round(float(line[0]),2)
				if pct > threshold_pct:
					species = line[-1]
					species = re.search('([A-Z].*)', species).group(1)
					to_report = species+':'+str(pct)
					species_conta.append(to_report)

	return(genus_conta, species_conta)


def Kraken2_contamination(f, expected_genus, expected_species, threshold_pct = 2):
	"""
	"""
	print("Expected genus: %s"%(expected_genus))
	print("Expected species: %s"%(expected_species))
	# get genus and species found in kraken2 output
	genus_conta, species_conta = kraken_extract(f, threshold_pct)

	# format results
	genus_level = {}
	for g in genus_conta:
		g = g.split(':')
		g_value = g[1]
		g = g[0]
		if g not in genus_level:
			genus_level[g]={'pct': g_value+'%', 'species': []}
	for s in species_conta:
		s = s.split(':')
		s_value = s[1]
		s = s[0]
		g = s.split(' ')[0]
		if g not in genus_level:
			genus_level[g]={'pct': s_value+'%', 'species': [s.split(' ')[1]+'('+s_value+'%)']}
		else:
			genus_level[g]['species'].append(s.split(' ')[1]+'('+s_value+'%)')

	to_report = []
	for g in genus_level:
		if g != expected_genus :
			to_report.append(g+' '+genus_level[g]['pct'])

		else:
			for s in genus_level[g]['species']:
				if not re.search('^'+expected_species, s) :
					to_report.append(g+' '+s)
	
	if to_report:
		return(' - '.join(to_report))
	else:
		return('False')


def main():

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
		
	Arguments=parser.parse_args()

	print('Start contamination analysis')
	print('')

	print('Search for confindr contamination')
	confindr_results = Confindr_contamination(Arguments.confindr,Arguments.confindr_SNVs_threshold)
	print('')

	print('Search for kraken2 contamination')
	kraken2_results = Kraken2_contamination(Arguments.kraken2, Arguments.genus, Arguments.species, threshold_pct = 2)
	print('')
	print('')


	print('Results:')
	with open(Arguments.output, 'w') as outf:
		outf.write('Confindr results: %s\n'%(confindr_results))
		print('Confindr results: %s'%(confindr_results))

		outf.write('Kraken2 results: %s'%(kraken2_results.strip()))
		print('Kraken2 results: %s\n'%(kraken2_results.strip()))

if __name__ == "__main__":
	main()		