import os, re

"""
Réplication des lignes de commandes linux:
head -n 1 {input} | sed 's/^FILE/#FILE/' > {output}
paste <(echo {wildcards.sample}) <(tail -n 1 {input} | cut -f 2- | sed 's/INF-//g' | sed -e 's/ALM\|ASM\|LNF\|NIPH\|NIPHEM\|PLOT[0-9]/-/g' | sed 's/NA/-/g' | sed -e '$a\\') >> {output}
"""

# Get sample id
sample_file = snakemake.input['assembly']
sample = os.path.basename(sample_file).replace('.fa', '')
print("sample name: %s"%(sample))

# Alleles results
allele_dir = os.path.dirname(snakemake.input['schema'])
allele = allele_dir+'/'+snakemake.params['genus']+'/results/results_alleles.tsv'




new, filtered = [], []
allele_line = []
with open(allele) as f:
    header = f.readline()
    for line in f:
        line = line.rstrip().split('\t')
        print(line[0])
        found = False
        if line[0] == sample:
            found = True
            name = line[0]
            allele_line.append(name)

            i=0
            for a in line[1:]:
                if re.search('\*', a) or re.search('INF', a):
                    new.append(i)
                    a = re.search('([0-9]*)$', a).group(1)
                elif re.search('LNF', a) or re.search('PLOT', a) or re.search('NIPH', a) or re.search('LOTSC', a) or re.search('ALM', a) or re.search('ASM', a):
                    a = '-'
                    filtered.append(i)
                else:
                    a = re.search('([0-9]*)$', a).group(1)
                allele_line.append(a)
                i+=1
            break
if not found:
    raise ValueError('A very bad thing happends!!! Sample is not found in chewbbaca output')
    
    
print('Number of new allele:', len(new))
print('Number of filtered allele:', len(filtered))

with open(snakemake.output['profile'], 'w') as out:
    out.write(header)
    out.write('\t'.join(allele_line)+'\n')
    
  
    
    