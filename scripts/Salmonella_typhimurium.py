#!python
import re
import pandas as pd

pattern = ">"
Mdh= "Mdh"
fliC = "fliC"
fliA_fliB = "fliA_fliB"
fljB = "fljB"



#Profils PCR attendus a comparer avec les profils obtenus

typhimirium1 = [ 'Mdh', 'fliC', 'fliA_fliB']
incoherent = [ 'Mdh', 'fliC', 'fljB','fliA_fliB']
autre_serovar1 = ['fliC', 'fliA_fliB']
typhimirium2 = [ 'Mdh', 'fljB', 'fliA_fliB']
autre_serovar2 = ['fljB','fliA_fliB']
immobile = ['Mdh', 'fliA_fliB']
autre_serovar_immobile = ['fliA_fliB']



#1er : Lecture du fichier de sortie isPCR 
pcr_result = open(snakemake.input['isPCR_result'])
#pcr_result = open("test/Salmonella_enterica/Salmonella/typhivar/Salmonella_enterica_IsPCR.txt")

ID_line = []
for line in pcr_result:
  if re.search(pattern, line):
    ID_line.append(line)



#2eme : Lecture du fichier de sortie blast de l'amplicon fliA_fliB
try:
    blast_result = pd.read_csv(snakemake.input['blast_result'], delimiter='\t', header=None)
    #blast_result = pd.read_csv("test/Salmonella_enterica/Salmonella/typhivar/Salmonella_enterica_fliA_fliB_blastn.txt", delimiter='\t', header=None)
    couv_blast = max(abs(blast_result[3]-blast_result[4])/185)
    ident_blast = max(blast_result[5])
except:
    """
    Dans le cas où aucun résultats dans l'output du blast
    """
    couv_blast = 0
    ident_blast = 0

#Etablissement du profil PCR 

profil_PCR = []


## A partir des resultats IsPCR
for i in range(len(ID_line)):
  if re.search(Mdh, ID_line[i]):
    profil_PCR.append(Mdh)
  if re.search(fliC, ID_line[i]):
    profil_PCR.append(fliC)
  if re.search(fljB, ID_line[i]):
    profil_PCR.append(fljB)

##A partir des resultats Blast
if float(couv_blast) > 0.99:
  if float(ident_blast) > 99:
    profil_PCR.append(fliA_fliB)


#Interpretation de profils 

if profil_PCR == typhimirium1:
  serotype = 'variant confirme typhimurium'
if profil_PCR == incoherent:
  serotype = 'variant incoherent'
if profil_PCR ==  autre_serovar1:
  serotype = 'autre serovar que typhimurium'
if profil_PCR == typhimirium2:
  serotype = 'variant confirme typhimurium'
if profil_PCR == autre_serovar2:
  serotype = 'autre serovar que typhimurium'
if profil_PCR == immobile:
  serotype = 'Variant immobile de typhimurium'
if profil_PCR == autre_serovar_immobile:
  serotype = 'Variant immobile autre serovar'
if profil_PCR != typhimirium1:
  if profil_PCR != incoherent:
    if profil_PCR !=  autre_serovar1:
      if profil_PCR != typhimirium2:
        if profil_PCR != autre_serovar2:
          if profil_PCR != immobile:
            if profil_PCR != autre_serovar_immobile:
              serotype = 'non interpretable'
with open(snakemake.output['typhymurium'], 'w') as out:
#with open("toto", 'w') as out:
  out.write(serotype+'\n')