import os
import csv



CC_dict = {}
with open(snakemake.input['CC']) as f:
    for line in f:
        line = line.rstrip().split('\t')
        CC_dict[line[0]] = line[-2] ## ST = CC

with open(snakemake.input['mlst_report']) as f:
    for line in f:
        line = line.rstrip().split('\t')
        st = line[2]
        
        if st in CC_dict:
            cc = CC_dict[st]
        else:
            cc = '-'

with open(snakemake.output['CC_report'], 'w') as out:
    out.write(cc)