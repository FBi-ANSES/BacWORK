# snakemake -s snakefile.py --use-conda --cores 16 --forceall

import glob
import os 

# Sérotype associé à chaque combinaison de gènes
dict_serotype = {
    ("Prfa", "lmo0737", "Prs"): 'serotype IIa',
    ("Prfa", "ORF2819", "Prs"): 'serotype IIb',
    ("lmo1118", "lmo0737", "Prfa", "Prs"): 'serotype IIc',
    ("Prfa", "ORF2819", "ORF2110", "Prs"): 'serotype IVb'
}


# Liste des gènes présents dans chaque fichier BLAST
liste_gene = []

with open(snakemake.input['blast'], 'r') as f_in :
    for line in f_in:
        # Séparer les lignes en une liste de mots, puis extraire le 1er mot de chaque ligne
        mots = line.split()
        nom_gene = mots[0]
        nom_gene = nom_gene.split('~~')[0]
        
        # Ajouter le nom de chaque gène dans la liste             
        liste_gene.append(nom_gene)
        
with open(snakemake.output['serotype'], "w") as f_out:
    serotype = False    
    for key, value in dict_serotype.items(): # Parcourt chaque tuple du dictionnaire 
        if set(key) == set(liste_gene) :
            f_out.write(f"{value}")
            serotype = True
    
    if not serotype:
        f_out.write(f"undetermined") 
        
            


