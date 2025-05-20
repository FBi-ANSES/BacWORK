import os, sys
import argparse
"""
Le filtre de couverture est réalisé sur l'information de couverture données par shovill lors de l'assemblage
Cette info est donnée dans le header du contigs

exemple:
>contig00001 len=1470795 cov=102.6 corr=0 origname=NODE_1_length_1470795_cov_102.647599_pilon sw=shovill-spades/1.1.0 date=20211216
"""


fasta_in = snakemake.input['contigs']
plasmid = f'{snakemake.input["plascope_dir"]}/{snakemake.params["sample"]}_plasmid.fasta'
fasta_out = snakemake.output['contigs']
report_out = snakemake.output['report']


minCov = float(snakemake.params['mincov'])
smp = snakemake.params['sample']


# read fasta plasmid if exist and create a list of predicted plasmid contigs 
header_plasmid = []
if os.path.exists(plasmid):
	with open(plasmid, 'r') as file:
		lines = file.readlines()
		file.close()
		for line in lines :
			if line[0] == '>':
				node_nb = '_'.join(line[1:].split('_')[0:2])+'_'
				header_plasmid.append(node_nb)


with open(report_out, 'w') as report_out:
	header_to_report = ['contig_id', 'shovill_id', 'length', 'coverage']
	report_out.write('\t'.join(header_to_report)+'\n')
	
	with open(fasta_out, 'w') as fasta_out:
		contig_n = 1
		with open(fasta_in) as f:
			skip = False
			for line in f:
				if line[0] == '>': # found header
					cov = float(line.split(' ')[2].replace('cov=',''))
					
					if cov < minCov:
						skip = True
					else:
						skip = False
						plasmid_statut = ''					
						for node in header_plasmid:
							if node in line :
								plasmid_statut = '_predicted_plasmid'
								continue								
						line = line.strip().split(' ')
						new_header = smp+'_'+str(contig_n)+plasmid_statut
						fasta_out.write('>'+new_header+'\n')
						contig_n += 1
						to_report = [new_header, 
										line[0].replace('>', ''), 
										line[1].replace('len=', ''),
										line[2].replace('cov=', '')]
						report_out.write('\t'.join(to_report)+'\n')
				else:
					if not skip :
						fasta_out.write(line.strip().upper()+'\n')

