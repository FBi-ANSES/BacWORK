import sys, os, re
import json

genus_species = {
            'Listeria': 'monocytogenes',
            'Staphylococcus': 'aureus',
            'Salmonella': 'enterica',
            'Clostridium': 'perfringens'
            }
            
            
run = sys.argv[1]
try:
    genus = os.path.basename(os.path.dirname(os.path.dirname(run)))
    genus_species[genus]
except:
    genus = os.path.basename(os.path.dirname(run))
    genus_species[genus]
# else:
    # raise ValueError('error in genus and species attribution')


config_ = {
	"WORKDIR" : "/scratch/v.chesnais/BacWork/Validation_Salmonella",
	"FASTQ_DIR" : run,
	"samples": []
}




for s in os.listdir(run):
    if re.search('_R1.fastq.gz', s):
        s_name = s.replace('_R1.fastq.gz', '')
        
        s_dict = {
            "SampleID": s_name,
            "Project": "",
            "Supplier": "",
            "Sequencing_center": "",
            "Phylogeny": {
                "Genus": genus,
                "Species": genus_species[genus]
                }
            }
        
        config_['samples'].append(s_dict)

with open('config/config.json', 'w', encoding='utf-8') as f:
    json.dump(config_, f, ensure_ascii=False, indent=4)