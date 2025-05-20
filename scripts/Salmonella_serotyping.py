import pandas as pd
import os, csv, json, re
from collections import Counter


## PARSING FUNCTIONS
def remove_optional(serovar):
    try:
        serovar = serovar.split('|')
    except:
        serovar = []
    out_serovar = []
    for s in serovar:
        tmp = s.split(':')
        if len(tmp) == 4:
            s = ':'.join(tmp[:-1])
        out_serovar.append(s)
    return(out_serovar)
    
def sistr_output(sistr, out_writer):
    with open(sistr) as f:
        sistr = json.load(f)[0]
        print(sistr.keys())
    out_writer.writerow(['serovar_cgmlst', sistr['serovar_cgmlst']])
    out_writer.writerow(['serovar', sistr['serovar']])
    out_writer.writerow(['serovar_antigen', sistr['serovar_antigen']])
    out_writer.writerow(['mash_serovar', sistr['mash_serovar']])
    out_writer.writerow(['mash_subspecies', sistr['mash_subspecies']])
    
    out = remove_optional(sistr['serovar_cgmlst'])
    out += remove_optional(sistr['serovar'])
    out += remove_optional(sistr['serovar_antigen'])
    out += remove_optional(sistr['mash_serovar'])
    subpecies = remove_optional(sistr['mash_subspecies'])[0]
    return(out, subpecies)

def seqsero_output(seqsero, out_writer):
    seqsero = pd.read_csv(seqsero, sep='\t').to_dict()
    out_writer.writerow(['Predicted serotype', seqsero['Predicted serotype'][0]])
    out_writer.writerow(['Predicted identification', seqsero['Predicted identification'][0]])
    
    out = remove_optional(seqsero['Predicted serotype'][0])
    subpecies = seqsero['Predicted identification'][0]
    try:
        subpecies = re.search(' ([a-zA-Z]*) \(subspecies', subpecies).group(1)
    except:
        subpecies = subpecies
    return(out, subpecies)
    
    
## MAIN FUNCTIONS
def define_serovar(f_sistr, f_seqsero, f_seqsero_reads, out_writer):
    out_writer.writerow(['SeqSero2 - results'])
    seqsero, seqsero_S = seqsero_output(f_seqsero, out_writer)
    print(f_seqsero)
    out_writer.writerow([' '])
    out_writer.writerow(['SeqSero2 - reads - results'])
    seqseroR, seqseroR_S = seqsero_output(f_seqsero_reads, out_writer)
    print(f_seqsero_reads)
    out_writer.writerow([' '])
    out_writer.writerow(['Sistr - results'])
    sistr, sistr_S = sistr_output(f_sistr, out_writer)
    print(f_sistr)
    out_writer.writerow([' '])
    
    # define subspecies
    if sistr_S:
        serovar_subspecies = sistr_S
    else:
        serovar_subspecies = seqsero_S
    
    # search frequent serovar
    serovar = Counter(seqsero+sistr)
    max_occurence = max(serovar.values())
    catch = []
    for s in serovar:
        if serovar[s] == max_occurence:
            catch.append(s)
    if len(catch) == 1:
        final_serovar = catch[0]
    else:
        final_serovar = 'NA'
    print(final_serovar)
    return(final_serovar, serovar_subspecies)

#
# with open('toto.txt', 'w') as out_report:
with open(snakemake.output['txt'], 'w') as out_report:
    out_writer = csv.writer(out_report, delimiter='\t', lineterminator="\n")
    final_serovar, final_subspecies = define_serovar(snakemake.input['sistr'], snakemake.input['seqsero'], snakemake.input['seqseroReads'], out_writer)
    # final_serovar, final_subspecies = define_serovar("test/Salmonella_enterica/Salmonella/sistr/Salmonella_enterica_alleles.json", "test/Salmonella_enterica/Salmonella/seqsero2/SeqSero_result.tsv", "test/Salmonella_enterica/Salmonella/seqsero2Reads/SeqSero_result.tsv", out_writer)

if final_serovar.lower() == "typhimurium":
    
    # with open("test/Salmonella_enterica/Salmonella/typhivar/Salmonella_enterica_typhimurium_variant.txt") as f:
    with open(snakemake.input['typhymurium']) as f:
        annot = f.readline().rstrip()
        final_serovar = final_serovar + '(' + annot + ')'

# with open("toto2.txt", 'w') as out_report: 
with open(snakemake.output['serotyping'], 'w') as out_report:

    out_writer = csv.writer(out_report, delimiter='\t', lineterminator="\n")
    out_writer.writerow([final_serovar])
    out_writer.writerow([final_subspecies])