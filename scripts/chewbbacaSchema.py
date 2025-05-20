import os, shutil, glob
import time



# List all available assemblies
print('Find samples')
assemblies = {}
input_f = snakemake.input['fasta']
for f in input_f:
    assemblies[os.path.basename(f).replace('.fa', '')] = f


# Search genus expected for every assembly: copy chewbbaca reference files and dispatch assemblies per genus
print('Search genus to use for every samples...')
database_genus = {}
cgMLST_SPECIES_LIST = ["Staphylococcus","Salmonella","Listeria","Clostridium","Vibrio","Taylorella"]

for sample in snakemake.params['config_json']:
    print(sample)
    genus = sample['Phylogeny']['Genus']
    if genus not in cgMLST_SPECIES_LIST:
        continue
    if genus not in database_genus:
        print(genus)
        db_path = os.path.join(
                        os.path.dirname(snakemake.output['schema']),
                        genus)
        
        # Copy database to avoid permission error
        print('... search %s database...'%(genus))
        if not os.path.isdir(os.path.join(db_path, 'schema')):
            print('... %s database not found: start copy database'%(genus))
            shutil.copytree(snakemake.params['schema'][genus], os.path.join(db_path, 'schema'))
        database_genus[genus] = {'schema':os.path.join(db_path, 'schema'),
                                 'genome_list': os.path.join(db_path, 'genome_list.txt'),
                                 'sample' :[]
                                }
    database_genus[genus]['sample'].append(sample['SampleID'])



# Launch chewbbaca for all genus
for genus in database_genus:
    genus = database_genus[genus]
    out_path = os.path.dirname(genus['genome_list'])
    # print(out_path)
    # Write the genome_list input file with all assemblies
    with open(genus['genome_list'], 'w') as w:
        for s in genus['sample']:
            w.write(assemblies[s]+'\n')
    
    # Chewbbaca and rename results directory
    if os.path.isdir(out_path+'/results'):
        all_results = glob.glob(out_path+'/result*')
        # print(all_results)
        for result in all_results:
            shutil.rmtree(result)
    
    time.sleep(30)
    print('Launch chewBBACA for %s !!'%(genus))
    cmd = "chewBBACA.py AlleleCall -i "+genus['genome_list']+" -g "+genus['schema']+" -o "+out_path+" --cpu "+str(snakemake.threads)+" --gl "+genus['schema']+"/listGenes.txt"
    os.system(cmd)
    time.sleep(30)
    out_dir = glob.glob(out_path+'/results_*/results_alleles.tsv')
    shutil.move(os.path.dirname(out_dir[0]), out_path+'/results')
    time.sleep(30)

with open(snakemake.output['schema'], 'w') as w:
    w.write('end')