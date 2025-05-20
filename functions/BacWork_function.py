import shutil
import json
import glob
import os

cgMLST_SPECIES_LIST = ["Staphylococcus","Salmonella","Listeria","Clostridium","Vibrio","Taylorella"]


def makedirs(WORKDIR):
	if not os.path.exists(WORKDIR):
		os.makedirs(WORKDIR)

def get_genus(sample_id):
	for element in config["samples"]:
		if element["SampleID"] == sample_id:
			return element["Phylogeny"]["Genus"]
			
def get_species(sample_id,kraken_report):
	for element in config["samples"]:
		if element["SampleID"] == sample_id:
			if element["Phylogeny"]["Species"] == "unknown" :
				return get_species_from_kraken2(kraken_report,sample_id)
			else:
				return element["Phylogeny"]["Species"]
				
def get_Project(sample_id):
	for element in config["samples"]:
		if element["SampleID"] == sample_id:
			return element["Project"]

def get_Supplier(sample_id):
	for element in config["samples"]:
		if element["SampleID"] == sample_id:
			return element["Supplier"]	

def get_Sequencing_center(sample_id):
	for element in config["samples"]:
		if element["SampleID"] == sample_id:
			return element["Sequencing_center"]			

def grep_accession_from_mash_screen_output(file_path):
	with open(file_path, 'r') as f:
		line = f.readline()
	accession = line.split('\t')[4]
	return accession
		
def get_species_from_kraken2(kraken_report,sample_id):
	genus = get_genus(sample_id)
	with open(kraken_report, 'r') as f:
		lines=f.readlines()
		for line in lines:
			if f"{genus} " in line :
				return line.rstrip().split(' ')[-1]

def get_database_dir(sample, plascope_dir,kraken_report):
	genus = get_genus(sample)
	species = get_species(sample,kraken_report)
	database_dir = f"{plascope_dir}/{genus}_{species}/"
	return database_dir
	
def get_database_name(sample, plascope_dir,kraken_report):
	genus = get_genus(sample)
	species = get_species(sample,kraken_report)
	database_name = f"all_{genus}_{species}_plascopedb"
	return database_name

def get_reads_preprocessing_files_list(config):
	reads_preprocessing_files = []
	for sample in config["samples"]:
		#fastQC html files
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/fastqc/raw/{sample['SampleID']}_R1_fastqc.html")
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/fastqc/raw/{sample['SampleID']}_R2_fastqc.html")
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/fastqc/clean/{sample['SampleID']}_clean_1_fastqc.html")
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/fastqc/clean/{sample['SampleID']}_clean_2_fastqc.html")
		#contamination
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/contamination/report.tsv")
		#bbmap file
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/bbmap/covstats.txt")
		#speciesfinder file
		reads_preprocessing_files.append(f"{WORKDIR}/{sample['SampleID']}/reads_preprocessing/speciesFinder/data.txt")
	return reads_preprocessing_files

def get_assembly_and_annotation_files_list(config):
	assembly_and_annotation_files = []
	for sample in config["samples"]:
		#mauvecm file
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/mauvecm/contigs.fa.fas")
		#quast report
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/quast/report.tsv")
		#bakta
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/bakta/{sample['SampleID']}.gbff")
		#busco
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/busco/short_summary.json")
		#resfinder file
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/resfinder/ResFinder_results.txt")
		#virulencefinder file
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/virulencefinder/results_tab.tsv")
		#abricate_vfdb file
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/abricate/vfdb.tsv")
		#vibrant dir
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/vibrant/VIBRANT_phages.fa")
		#plasmidfinder file
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/plasmidfinder/data.json")
		#plascope dir
		assembly_and_annotation_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/plascope/{sample['SampleID']}_PlaScope/PlaScope_predictions")

	return assembly_and_annotation_files

def get_assembly_list(config):
	assembly_files = []
	for sample in config["samples"]:
		if get_genus(sample['SampleID']) in cgMLST_SPECIES_LIST:
			assembly_files.append(f"{WORKDIR}/{sample['SampleID']}/assembly_and_annotation/filter/{sample['SampleID']}.fa")
	return assembly_files 
	
	
def get_sequence_typing_files_list(config):
	sequence_typing_files = []
	for sample in config["samples"]:
		if get_genus(sample['SampleID']) in cgMLST_SPECIES_LIST:
			sequence_typing_files.append(f"{WORKDIR}/{sample['SampleID']}/sequence_typing/mlst/{sample['SampleID']}_mlst.txt")
			sequence_typing_files.append(f"{WORKDIR}/{sample['SampleID']}/sequence_typing/chewbbaca/alleles_hash.tsv")
			sequence_typing_files.append(f"{WORKDIR}/{sample['SampleID']}/sequence_typing/chewbbaca/alleles.tsv")
	return sequence_typing_files 	

	
def get_species_specific_files_list(config):	
	species_specific_files = []
	for sample in config["samples"]:
		if get_genus(sample['SampleID'])=="Bacillus":
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Bacillus/btyper3/{sample['SampleID']}_final_results.txt")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Bacillus/Bt_detect/{sample['SampleID']}.tsv")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Bacillus/abricate/VBC.tsv")
		if get_genus(sample['SampleID'])=="Escherichia":
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Escherichia/clermontyping/{sample['SampleID']}_phylogroups.txt")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Escherichia/stecfinder/{sample['SampleID']}_stecfinder.txt")
		if get_genus(sample['SampleID'])=="Campylobacter":	
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Campylobacter/abricate/host_markers.tsv")
		if get_genus(sample['SampleID'])=="Salmonella":		
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Salmonella/seqsero2/SeqSero_result.txt")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Salmonella/serotyping/{sample['SampleID']}.txt")
		if get_genus(sample['SampleID'])=="Staphylococcus":
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Staphylococcus/spatyper/{sample['SampleID']}_spatyper.txt")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Staphylococcus/staphopia-sccmec/{sample['SampleID']}_staphopia-sccmec.txt")
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Staphylococcus/naura/matrix.tsv")
		if get_genus(sample['SampleID'])=="Listeria":
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Listeria/serotyping/{sample['SampleID']}_serotyping.txt")
		if get_genus(sample['SampleID'])=="Clostridium":	
			species_specific_files.append(f"{WORKDIR}/{sample['SampleID']}/Clostridium/abricate/virulence_markers.tsv")	
	return species_specific_files 		

def multiqc_input_files(config):
	reads_preprocessing_files = get_reads_preprocessing_files_list(config)
	assembly_and_annotation_files = get_assembly_and_annotation_files_list(config)
	return reads_preprocessing_files + assembly_and_annotation_files
	
def clean_tmp_files():
	for folder in glob.glob(f"{WORKDIR}/cgMLST/*/schema"):
		if os.path.isdir(folder):
			shutil.rmtree(folder)
	for folder in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/busco/analyse_busco"):
		if os.path.isdir(folder):
			shutil.rmtree(folder)
	for folder in glob.glob(f"{WORKDIR}/*/Escherichia/clermontyping/analysis_*"):
		if os.path.isdir(folder):
			shutil.rmtree(folder)	
	for folder in glob.glob(f"{WORKDIR}/*/Escherichia/clermontyping/ClermonTyping"):
		if os.path.isdir(folder):
			shutil.rmtree(folder)			
	for file in glob.glob(f"{WORKDIR}/*/reads_preprocessing/kraken2/out.krepport"):
		if os.path.isfile(file):
			os.remove(file)
	for file in glob.glob(f"{WORKDIR}/*/reads_preprocessing/mash_screen/reference.fa.sslist"):
		if os.path.isfile(file):
			os.remove(file)
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/shovill/spades.fasta"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/shovill/contigs.gfa"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/shovill/shovill.log"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/mauvecm/contigs.fa.fas.sslist"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/mauvecm/reference.fa.sslist"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/mauvecm/reference.fa"):
		if os.path.isfile(file):
			os.remove(file)			
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/quast/report.html"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/bakta/{{sample}}.embl"):
		if os.path.isfile(file):
			os.remove(file)	
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/bakta/{{sample}}.json"):
		if os.path.isfile(file):
			os.remove(file)		
	for file in glob.glob(f"{WORKDIR}/*/assembly_and_annotation/busco/{{sample}}.fa"):
		if os.path.isfile(file):
			os.remove(file)		
	for file in glob.glob(f"{WORKDIR}/*/Staphylococcus/naura/bakta.gbk"):
		if os.path.isfile(file):
			os.remove(file)					
			
def resume_files_list(config):
	resume_files = []
	for sample in config["samples"]:
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/summary.json")
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/tools_list.tsv")
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/db_version.tsv")
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/quality_threshold.json")
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/config_tools.json")
		resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/{sample['SampleID']}_assembly.fa")
		if get_genus(sample['SampleID'])=="Salmonella":	
			resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/resume_Salmonella")
		if get_genus(sample['SampleID'])=="Bacillus":	
			resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/resume_Bacillus")
		if get_genus(sample['SampleID'])=="Staphylococcus":	
			resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/resume_Staphylococcus")
		if get_genus(sample['SampleID'])=="Escherichia":	
			resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/resume_Escherichia")
		if get_genus(sample['SampleID'])=="Listeria":	
			resume_files.append(f"{WORKDIR}/{sample['SampleID']}/resume/resume_Listeria")		
	return resume_files 	

def get_all_input(config):
	reads_preprocessing_files = get_reads_preprocessing_files_list(config)
	assembly_and_annotation_files = get_assembly_and_annotation_files_list(config)
	sequence_typing_files = get_sequence_typing_files_list(config)
	species_specific_files = get_species_specific_files_list(config)
	resume_files = resume_files_list(config)
	final_input_list = []
	
	if config["Analysis"]["reads preprocessing"] == "yes":
		final_input_list = final_input_list + reads_preprocessing_files 
	if config["Analysis"]["assembly_and_annotation"] == "yes":
		final_input_list = final_input_list + assembly_and_annotation_files
	if config["Analysis"]["sequence_typing"] == "yes":
		final_input_list = final_input_list + sequence_typing_files	
	if config["Analysis"]["species_specific_steps"] == "yes":
		final_input_list = final_input_list + species_specific_files
	if config["Analysis"]["resume"] == "yes":		
		final_input_list = final_input_list + resume_files 
		final_input_list.append(f"{WORKDIR}/Bacwork_reports/multiqc_report.html")
		final_input_list.append(f"{WORKDIR}/Bacwork_reports/resume.json")
		final_input_list.append(f"{WORKDIR}/Bacwork_reports/resume.tsv")
	return final_input_list
