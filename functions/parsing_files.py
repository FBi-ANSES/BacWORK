import argparse
import os, sys
import json
import datetime
import uuid
import glob
from shutil import copyfile

####################################################################################
#  FASTQ HEADER FILE  ##############################################################
def detect_SequencingTechnology(reads_filepath):

    tmpFile = str(uuid.uuid1())
    cmd = f"zcat {reads_filepath} | head -n 1 > {tmpFile}"
    os.system(cmd)

    with open(tmpFile, "r") as f:
        for line in f.readlines():
        
            if line[0:6] == "@HWI-M" or line[0:2] == "@M":
                techno = 'MiSeq'

            elif line[0:6] == "@HWI-C" or line[0:2] == "@C":
                techno = 'HiSeq 1500'
                
            elif line[0:6] == "@HWUSI" :
                techno = 'GAIIx'

            elif line[0:6] == '@HWI-D' or line[0:2]  == '@D':
                techno = 'HiSeq 2500'
                
            elif line[0:2]  == '@J':
                techno = 'HiSeq 3000'
                
            elif line[0:2]  == '@K':
                techno = 'HiSeq 3000/4000'  
                
            elif line[0:2]  == '@E':
                techno = 'HiSeq X'  
                
            elif line[0:2] == '@N' :
                techno = 'NextSeq'
                
            elif line[0:2] == '@A' :
                techno = 'NovaSeq'
                
            elif line[0:3] == '@MN' :
                techno = 'MiniSeq'
                
            else :
                techno = "illumina"
            break
    
    os.remove(tmpFile) 
    return techno


####################################################################################
#  FROM CONTAMINATION REPORT  ######################################################
def get_kraken2_status(contamination_report):
    with open(contamination_report, "r") as f:
        conta = f.readlines()
        kraken2 = conta[1].split(' ')[2].rstrip()
    return(kraken2)

####################################################################################
#  FROM FASTP JSON  ################################################################
def get_number_of_reads_before_filtering(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file) 
        return str(data["summary"]["before_filtering"]["total_reads"])

def get_reads_len_before_filtering(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file) 
        return str(data["summary"]["before_filtering"]["read1_mean_length"])

def get_number_of_reads_after_filtering(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file) 
        return str(data["summary"]["after_filtering"]["total_reads"])

def get_reads_len_after_filtering(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file) 
        return str(data["summary"]["after_filtering"]["read1_mean_length"])

def get_q30_before_filtering(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file) 
        return str(data["summary"]["before_filtering"]["q30_rate"])

def get_percent_of_filtered_reads(fastp_json):
    with open(fastp_json) as json_file:
        data = json.load(json_file)
        _ = data["summary"]["after_filtering"]["total_reads"] / data["summary"]["before_filtering"]["total_reads"]
        return str(round(_*100, 2))


####################################################################################
#  FROM BBMAP REPORT  ##############################################################
def get_closer_ref_id(covstats_file):
    with open(covstats_file, "r") as f:
        line = f.readlines()[1].rstrip().split('\t')
        return line[0].split(' ')[0]

def get_DeepCoverage_closerRef(covstats_file):
    with open(covstats_file, "r") as f:
        line = f.readlines()[1].strip().split('\t')
        return line[1]

def get_DeepCoverageStd_closerRef(covstats_file):
    with open(covstats_file, "r") as f:
        line = f.readlines()[1].strip().split('\t')
        return line[10]

def get_BreadthCoverage_closerRef(covstats_file):   
    with open(covstats_file, "r") as f:
        line = f.readlines()[1].strip().split('\t')
        return line[4]

def get_closer_ref_length(covstats_file):
    with open(covstats_file, "r") as f:
        line = f.readlines()[1].strip().split('\t')
        return line[2]


####################################################################################
#  FROM QUAST REPORT  ##############################################################
def get_NbContigs(quast_report):
    with open(quast_report, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            if line[0] == "# contigs (>= 0 bp)":
                return line[1]

def get_TotalLength(quast_report):
    with open(quast_report, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            if line[0] == "Total length (>= 0 bp)":
                return line[1]

def get_LargestContig(quast_report):
    with open(quast_report, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            if line[0] == "Largest contig":
                return line[1]

def get_N50(quast_report):
    with open(quast_report, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            if line[0] == "N50":
                return line[1]

def get_GenomeFraction(quast_report):
    with open(quast_report, "r") as f:
        for line in f.readlines():
            line = line.rstrip().split('\t')
            if line[0] == "Genome fraction (%)":
                return line[1]


####################################################################################
#  FROM MLST   #####################################################################
def get_ST(mlst_report):
    with open(mlst_report, "r") as f:
        return f.readlines()[0].strip().split('\t')[2]


####################################################################################
#  FROM ASSEMBLY   #################################################################
def get_plasmid_contigs(assembly):
    nb_plasmid_contigs = 0
    with open(assembly, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "_predicted_plasmid" in line :
                nb_plasmid_contigs += 1
    return str(nb_plasmid_contigs)


####################################################################################
#  QUALITY TESTS   #################################################################
def quality_Reads_deepcoverage(deep_coverage):
    if float(deep_coverage) >= float(config["min_reads_deep_coverage"]):
        return "PASS"
    else:
        return "FALL"
    
def quality_Reads_breadthcoverage(breadth_coverage):
    if float(breadth_coverage) >= float(config["min_reads_breadth_coverage"]):
        return "PASS"
    else:
        return "FALL"
        
def quality_nbcontigs(nb_contigs):
    if float(nb_contigs) <= float(config["max_number_of_contigs"]):
        return "PASS"
    else:
        return "FALL"       
    
def quality_assembly_breadthcoverage(breadth_coverage):
    if float(breadth_coverage) >= float(config["min_assembly_breadth_coverage"]):
        return "PASS"
    else:
        return "FALL"

def get_conta_statut(contamination_report):
    with open(contamination_report, "r") as f:
        conta = f.readlines()
        confindr = conta[0].split(' ')[2].rstrip()
        kraken2 = conta[1].split(' ')[2].rstrip()
        if (confindr == "True") or (kraken2 != "False"):
            return "FALL"
        elif (confindr == "False") and (kraken2 == "False"):
            return "PASS"
        else:
            return "ERROR"


####################################################################################
#  SPECIES SPECIFIC OUTPUT   #######################################################
def get_serotype_from_seqsero(seqsero):
    with open(seqsero, "r") as f:
        line = f.readlines()[1].split('\t')
        antigenic_profile = line[7]
        serovar = line[8]
        return antigenic_profile, serovar

def get_serotype_from_sistr(sistr):
    with open(sistr, "r") as f:
        lines = f.readlines()
        match = re.search(r'"serovar":\s*"([^"]+)"', str(lines))
        return  match.group(1)
		
def get_bacillus_cluster(bt_detect):
    with open(bt_detect, "r") as f:
        return f.readlines()[0].rstrip().split('\t')[-1]

def get_btyper3_group(btyper):  
    with open(btyper, "r") as f:
        line = f.readlines()[1].rstrip().split('\t')
        group = line[17]
        taxon = line[18]
        return group, taxon
        
def get_spatype(spatyper):
    with open(spatyper, "r") as f:
        return f.readlines()[1].rstrip().split('\t')[-1]

def get_Escherichia_serotype(stecfinder):
    with open(stecfinder, "r") as f:
        line = f.readlines()[1].split('\t')
        serotype = line[3]
        if "Non-STEC" in line[9]:
            stec = "No"
        else :
            stec = "Yes"
        return serotype, stec

def get_Staphylococcus_spatype(spatyper):
    with open(spatyper, "r") as f:
        return f.readlines()[1].rstrip().split('\t')[-1]

def get_Listeria_serotype(serotype):
    with open(serotype, "r") as f:
        return f.readlines()[0].rstrip()
           
def get_Listeria_CC(CC_report):
    with open(CC_report, "r") as f:
        return f.readlines()[0].strip()
        


####################################################################################
#  CREATE RESUME OUTPUT   ##########################################################     
def make_resume_json(sample_id,output_file,fastp_json,clean_R1,bbmap_covstat,\
    quast_report,kraken_report,mlst_report,assembly,contamination_report):

    json_dico = {
        "BacWork_version":get_bacWork_version(),
        "SampleID":sample_id,
        "Supplier":get_Supplier(sample_id),
        "Project":get_Project(sample_id),
        "ProcessingDate":str(datetime.date.today()),
        "Phylogeny":{
            "Genus":get_genus(sample_id),
            "Species":get_species(sample_id,kraken_report),
            "ST":get_ST(mlst_report),
        },
        "Reads":{
            "Center":get_Sequencing_center(sample_id),
            "Predicted_Technology":detect_SequencingTechnology(clean_R1),
            "Reads_before_filtering":get_number_of_reads_before_filtering(fastp_json),
            "Reads_length_before_filtering":get_reads_len_before_filtering(fastp_json),
            "Reads_after_filtering":get_number_of_reads_after_filtering(fastp_json),
            "Q30_before_filtering":get_q30_before_filtering(fastp_json),
            "Reads_percent_filtered":get_percent_of_filtered_reads(fastp_json),
            "Reads_length_after_filtering":get_reads_len_after_filtering(fastp_json),
            "CloserRef":get_closer_ref_id(bbmap_covstat),
            "CloserRef_length":get_closer_ref_length(bbmap_covstat),
            "DeepCoverage_closerRef":get_DeepCoverage_closerRef(bbmap_covstat),
            "DeepCoverageStd_closerRef":get_DeepCoverageStd_closerRef(bbmap_covstat),
            "BreadthCoverage_closerRef":get_BreadthCoverage_closerRef(bbmap_covstat),
            "Kraken2_output":get_kraken2_status(contamination_report)
        },
        "Assembly":{
            "NbContigs":get_NbContigs(quast_report),
            "TotalLength":get_TotalLength(quast_report),
            "LargestContig":get_LargestContig(quast_report),
            "N50":get_N50(quast_report),
            "GenomeFraction":get_GenomeFraction(quast_report),
            "Plasmid_contigs":get_plasmid_contigs(assembly)
        },
    }
    
    json_dico["Quality_filter"] = {
        "Reads_deepcoverage":quality_Reads_deepcoverage(json_dico["Reads"]["DeepCoverage_closerRef"]),
        "Reads_breadthcoverage":quality_Reads_breadthcoverage(json_dico["Reads"]["BreadthCoverage_closerRef"]),
        "Contamination":get_conta_statut(contamination_report),
        "Nb_contigs":quality_nbcontigs(json_dico["Assembly"]["NbContigs"]),
        "Assembly_breadth_coverage":quality_assembly_breadthcoverage(json_dico["Assembly"]["GenomeFraction"])
    }
        
    with open(output_file, "w") as outfile:  
        json.dump(json_dico, outfile, indent = 4)
        
def add_Salmo_informations(resume_json, seqsero, sistr):    

    with open(resume_json) as json_file:
        data = json.load(json_file) 
        Seqsero2_antigenic_profile, Seqsero2_serovar = get_serotype_from_seqsero(seqsero)
        data["Phylogeny"]["Seqsero2_antigenic_profile"] = Seqsero2_antigenic_profile
        data["Phylogeny"]["Seqsero2_serovar"] = Seqsero2_serovar
        sistr_serovar = get_serotype_from_sistr(sistr)
        data["Phylogeny"]["Sistr_serovar"] = sistr_serovar
        
    with open(resume_json, "w") as outfile:  
        json.dump(data, outfile, indent = 4)    
        
def add_Bacillus_informations(resume_json, bt_detect, btyper):  

    with open(resume_json) as json_file:
        data = json.load(json_file) 
        group, taxon = get_btyper3_group(btyper)
        data["Phylogeny"]["group"] = group
        data["Phylogeny"]["taxon"] = taxon
        data["Phylogeny"]["bt_cluster"] = get_bacillus_cluster(bt_detect)
        
    with open(resume_json, "w") as outfile:  
        json.dump(data, outfile, indent = 4)    
        
def add_Escherichia_informations(resume_json, stecfinder):

    with open(resume_json) as json_file:
        data = json.load(json_file) 
        serotype, stec = get_Escherichia_serotype(stecfinder)
        data["Phylogeny"]["Serotype"] = serotype
        data["Phylogeny"]["STEC"] = stec
        
    with open(resume_json, "w") as outfile:  
        json.dump(data, outfile, indent = 4)

def add_Staph_informations(resume_json, spatyper):

    with open(resume_json) as json_file:
        data = json.load(json_file) 
        spatype = get_Staphylococcus_spatype(spatyper)
        data["Phylogeny"]["spatype"] = spatype
        
    with open(resume_json, "w") as outfile:  
        json.dump(data, outfile, indent = 4)
        
def add_Listeria_informations(resume_json, serotype, CC_report):

    with open(resume_json) as json_file:
        data = json.load(json_file) 
        data["Phylogeny"]["Serotype"] = get_Listeria_serotype(serotype)
        data["Phylogeny"]["CC"] = get_Listeria_CC(CC_report)
        
    with open(resume_json, "w") as outfile:  
        json.dump(data, outfile, indent = 4)        
            
def merge_summary_json(json_list, output_file):
    dico = {}
    for element in json_list:
        if "summary.json" in element:
            with open(element) as json_file:
                data = json.load(json_file)
                dico[data["SampleID"]] = data       
    with open(output_file, "w") as outfile: 
        json.dump(dico, outfile, indent = 4)


def json_to_tab(json_path, output_file):

    out = open(output_file,'w')
    
    tab_header_element = ["SampleID","Supplier","Project",\
    "Phylogeny.Genus","Phylogeny.Species","Phylogeny.ST",\
    # Salmonella specific  \
    "Phylogeny.Se_Seqsero2_antigenic_profile","Phylogeny.Se_Seqsero2_serovar",\ 
    "Phylogeny.Se_Sistr_serovar",\  
    # Bacillus specific \
    "Phylogeny.Bc_group","Phylogeny.Bc_taxon","Phylogeny.Bt_cluster",\
    # Escherichia specific \
    "Phylogeny.Ec_Serotype","Phylogeny.Ec_STEC",\
    # Listeria specific \
    "Phylogeny.Lm_Serotype",\
    "Phylogeny.Lm_CC",\
    # Staphylococcus specific \
    "Phylogeny.Sa_spatype",\
    "Reads.center","Reads.Predicted_Technology","Reads.Reads_before_filtering",\
        "Reads.Reads_length_before_filtering","Reads.Reads_after_filtering",\
        "Reads.Reads_length_after_filtering","Reads.CloserRef",\
        "Reads.DeepCoverage_closerRef","Reads.DeepCoverageStd_closerRef",\
        "Reads.BreadthCoverage_closerRef",\ 
    "Assembly.NbContigs","Assembly.TotalLength","Assembly.LargestContig",\
        "Assembly.N50","Assembly.GenomeFraction","Assembly.Plasmid_contigs",\
    "Quality.Reads_deepcoverage","Quality.Reads_breadthcoverage",\
    "Quality.Contamination","Quality.Nb_contigs","Quality.Assembly_breadth_coverage",\
    "ProcessingDate"]
    
    out.write('\t'.join(tab_header_element))
    
    with open(json_path) as json_file:
        data = json.load(json_file)
        header_list_element = []
        for sample in data :
            to_write = [data[sample]["SampleID"],\
                        data[sample]["Supplier"],\
                        data[sample]["Project"],\
                        data[sample]["Phylogeny"]["Genus"],\
                        data[sample]["Phylogeny"]["Species"],\
                        data[sample]["Phylogeny"]["ST"]\
                        ]
            out.write('\n'+'\t'.join(to_write))
            
            ### Salmonella
            if data[sample]["Phylogeny"]["Genus"] == "Salmonella":
                Se_specific = [data[sample]["Phylogeny"]["Seqsero2_antigenic_profile"],\
                                data[sample]["Phylogeny"]["Seqsero2_serovar"],\
                                data[sample]["Phylogeny"]["Sistr_serovar"]]
            else:
                Se_specific = ["-","-","-"]
            out.write('\t'+'\t'.join(Se_specific))

            ### Bacillus
            if data[sample]["Phylogeny"]["Genus"] == "Bacillus":
                Bc_specific = [data[sample]["Phylogeny"]["group"],\
                                data[sample]["Phylogeny"]["taxon"],\
                                data[sample]["Phylogeny"]["bt_cluster"]]
            else:
                Bc_specific = ["-","-","-"]
            out.write('\t'+'\t'.join(Bc_specific))

            ### Escherichia
            if data[sample]["Phylogeny"]["Genus"] == "Escherichia":
                Ec_specific = [data[sample]["Phylogeny"]["Serotype"],\
                                data[sample]["Phylogeny"]["STEC"]]
            else:
                Ec_specific = ["-","-"]
            out.write('\t'+'\t'.join(Ec_specific))
            
            ### Listeria
            if data[sample]["Phylogeny"]["Genus"] == "Listeria":
                Lm_specific = [data[sample]["Phylogeny"]["Serotype"],
                                data[sample]["Phylogeny"]["CC"]]
            else:
                Lm_specific = ["-", "-"]
            out.write('\t'+'\t'.join(Lm_specific))
            
            ### Staphylococcus
            if data[sample]["Phylogeny"]["Genus"] == "Staphylococcus":
                Sa_specific = [data[sample]["Phylogeny"]["spatype"]]
            else:
                Sa_specific = ["-"]
            out.write('\t'+'\t'.join(Sa_specific))

            to_write = [data[sample]["Reads"]["Center"],\
                        data[sample]["Reads"]["Predicted_Technology"],\
                        data[sample]["Reads"]["Reads_before_filtering"],\
                        data[sample]["Reads"]["Reads_length_before_filtering"],\
                        data[sample]["Reads"]["Reads_after_filtering"],\
                        data[sample]["Reads"]["Reads_length_after_filtering"],\
                        data[sample]["Reads"]["CloserRef"],\
                        data[sample]["Reads"]["DeepCoverage_closerRef"],\
                        data[sample]["Reads"]["DeepCoverageStd_closerRef"],\
                        data[sample]["Reads"]["BreadthCoverage_closerRef"]]
            out.write('\t'+'\t'.join(to_write))     

            to_write = [data[sample]["Assembly"]["NbContigs"],\
                        data[sample]["Assembly"]["TotalLength"],\
                        data[sample]["Assembly"]["LargestContig"],\
                        data[sample]["Assembly"]["N50"],\
                        data[sample]["Assembly"]["GenomeFraction"],\
                        data[sample]["Assembly"]["Plasmid_contigs"]]
            out.write('\t'+'\t'.join(to_write))

            to_write = [data[sample]["Quality_filter"]["Reads_deepcoverage"],\
                        data[sample]["Quality_filter"]["Reads_breadthcoverage"],\
                        data[sample]["Quality_filter"]["Contamination"],\
                        data[sample]["Quality_filter"]["Nb_contigs"],\
                        data[sample]["Quality_filter"]["Assembly_breadth_coverage"]]
            out.write('\t'+'\t'.join(to_write))
            
            out.write('\t'+data[sample]["ProcessingDate"])
             
    out.close()


####################################################################################
#  VERSION TRACKING   ##############################################################              
def make_tools_list(outfile, resfinder_version_file):

    out = open(outfile,'w')
    yaml_files = glob.glob("envs/*.yaml")
    for file in yaml_files :
        env = file.split('/')[-1].split('.')[0]
        out.write(f"############## {env} ##############\n")
        f = open(file,'r')
        lines = f.readlines()
        f.close()
        for line in lines :
            if '=' not in line :
                continue
            line = line.rstrip().split(' ')[-1]
            tool = line.split('=')[0]
            version = line.split('=')[-1]
            out.write(f"{tool}\t{version}\n")
        if "spatyper" in env :
            out.write(f"spatyper\t0.3.3\n")
        out.write(f"\n")    
        
    out.write(f"############## ResFinder ##############\n")
    res_file = open(resfinder_version_file,'r')
    res_version = res_file.readlines()[0].rstrip()
    res_file.close()
    out.write(f"ResFinder\t{res_version}\n")
    
    out.close()

'''
# DEPRECIATED
def make_db_version_file(outfile):  
    
    out = open(outfile,'w')
    yaml = open("config/config_path.json")
    lines = yaml.readlines()
    yaml.close()
    flag = False
    for line in lines :
        if "cgmlst_database" in line :
            flag = True
            continue
        if '/' in line :
            db = line.split('"')[1]
            version = line.split('/')[-2] # il doit y avoir un '/' à la fin de chaque path dans le fichier!
            if flag:
                out.write("cg/wgMLST database:")
            out.write(f"{db}\t{version}\n")
    out.close()     
'''            

def make_db_version_file(outfile):  
    
    out = open(outfile,'w')
    yaml = open("config/config_path.json")
    lines = yaml.readlines()
    yaml.close()
    flag = False
    for line in lines :
        if '/' in line :
            db = line.split('"')[1]
            version = "unkown"
            db_path = line.split('"')[3]
            if db == 'complete_genome_sketch_database' :
                db_path = '/'.join(db_path.split('/')[0:-1])
            if os.path.isfile(f"{db_path}/VERSION"):
                file = open(f"{db_path}/VERSION",'r')
                version = file.readlines()[0].rsplit()[0]
            out.write(f"{db}\t{version}\n")
    out.close()     		
		
			
def get_bacWork_version():
    try:
        version = os.popen("git describe --exact-match --abbrev=0").readlines()[0].replace("\n","")
    except:
        version = "unkown"
    return version  
    
    