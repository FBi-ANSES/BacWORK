#!/global/conda/envs/artwork/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys
import argparse

def get_parser():
	"""
	Parse arguments
	@return: arguments list
	@rtype: parser object
	"""
	parser = argparse.ArgumentParser(description='parse Fluidigm output and assign serovar')

	parser.add_argument('-i', action="store", dest='genome_list',
						nargs='+', required=True, help='genome assemblies (REQUIRED)')

	parser.add_argument('-t', action="store", dest='decision_table',
						type=str, default="table_bt.txt", help='Decision table (default:table_bt.txt)')
						
	parser.add_argument('-min_id', action="store", dest='min_identity',
						type=str, default='90', help='Delta max value for filter (default:90)')					

	parser.add_argument('-min_cov', action="store", dest='min_coverage',
						type=str, default="90", help='Output file (default:90)')	
						
	parser.add_argument('-db_dir', action="store", dest='db_dir',
						type=str, default="abricate_db", help='database dir (default:abricate_db)')	

	parser.add_argument('-T', action="store", dest='threads',
						type=str, default="4", help='number of threads (default:4)')								
						
	return parser


		
		
def decision_table_parser(file) :	

	f = open(file, 'r')
	lines = f.readlines()
	f.close()
	
	flag = False
	dico_marker = {}
	header = True
	header_list = []
	dico_association = {}
	result_key = {}
	i = 0
	
	for line in lines :
		if '###' in line:
			flag = True
			continue
		if not flag :
			line = line.rstrip().split('\t')
			dico_marker[line[0]] = line[1:]
		if flag and header :
			header = False
			header_list = line.rstrip().split('\t')[1:]
		elif flag :
			line = line.rstrip().split('\t')
			dico_association[i] = line[1:]
			if line[0] not in result_key :
				result_key[line[0]] = []
			result_key[line[0]].append(i)
			i+=1
			
	return 	dico_marker, header_list, dico_association,result_key	
		
		
def abricate_parser(abricate_file, mincov, minid):
	
	f = open(abricate_file, 'r')
	lines = f.readlines()
	f.close()

	marker_list = []
	header = True

	for line in lines :
		if header :
			header = False
			continue
		line = line.rstrip().split('\t')
		
		if line[5] not in marker_list and float(line[9])>=float(mincov) and float(line[10])>=float(minid):
			marker_list.append(line[5])
			
	return marker_list
		
		
def seqLen(fasta):

	file = open(fasta,'r')
	lines = file.readlines()
	file.close()
	
	dico_sequence = {}
	id = ''
	seqlen = 0
	for line in lines :
		line = line.rstrip()
		if len(line) == 0 :
			continue
		if line[0]=='>' :
			if id != '' :
				dico_sequence[id] = seqlen
			id = line[1:]
			seqlen = 0
		else :
			seqlen = seqlen + len(line)
	dico_sequence[id] = seqlen
	
	return dico_sequence
	
	
def blast_parser(blast_file, mincov, minid, dico_len):

	file = open(blast_file,'r')
	lines = file.readlines()
	file.close()
	
	marker_list = []
	
	for line in lines :
		line = line.rstrip().split('\t')
		if float(line[4]) < float(minid):
			continue
		cov = float(line[3])/float(dico_len[line[1]])*100
		if float(cov) < float(mincov):
			continue
		if line[1] not in marker_list :
			marker_list.append(line[1])
	
	return marker_list
	
	
	
def assign(marker_found,dico_marker,header_list,dico_association,genome_id,result_key):
	
	group_list = []
	result = genome_id
	for group in dico_marker :
		for marker in marker_found :
			if marker in dico_marker[group]:
				group_list.append(group)
				result = f'{result}\t{group}(+)'
				break
		if (group + '(') not in result :
			result = f'{result}\t{group}(-)'
					

	for Bacillus in dico_association :
		i = 0
		flag = True
		for element in dico_association[Bacillus]:
			if element == '1' and header_list[i] not in group_list :
				flag = False
				break
			elif element == '-1' and header_list[i] in group_list :
				flag = False
				break
			i+=1
		if flag :
			for element in result_key :
				if Bacillus in result_key[element]:
					result = f'{result}\t{element}'
					break
	
	return result
	
	
#main function	
def main():	

	
	# Get arguments 
	parser=get_parser()
	
	# Print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	Arguments=parser.parse_args()	
	
	# Process decision table
	dico_marker, header_list, dico_association,result_key = decision_table_parser(Arguments.decision_table)
	dico_len = seqLen(Arguments.db_dir)
	
	for assembly in Arguments.genome_list :
		genome_id = os.path.splitext(os.path.basename(assembly))[0]
		#abricate_out = genome_id + "_abricate.tsv"
		blast_out = genome_id + "_blast.tsv"
		#cmd = f"abricate --db {Arguments.db_name} --datadir {Arguments.db_dir} --threads {Arguments.threads} --quiet --minid {Arguments.min_identity} --mincov {Arguments.min_coverage} {assembly} > {abricate_out}"
		#print(cmd)
		cmd = f"blastx -query {assembly} -db {Arguments.db_dir} -out {blast_out} -outfmt '6 qseqid sseqid evalue length pident' -num_threads {Arguments.threads}"
		os.system(cmd)
		#marker_list = abricate_parser(abricate_out,Arguments.min_coverage,Arguments.min_identity)
		marker_list = blast_parser(blast_out,Arguments.min_coverage,Arguments.min_identity, dico_len)
		result = assign(marker_list,dico_marker,header_list,dico_association,genome_id,result_key)
		print(result)
		
		
if __name__ == "__main__":
	main()	   	
	