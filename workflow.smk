# Chargement des fichiers de configuration
configfile: "config/config.json"
configfile: "config/config_tools.json"
configfile: "config/config_path.json"
configfile: "config/config_quality.json"

# Chargement des fonction et variables globales
WORKDIR = config["WORKDIR"]
FASTQ_DIR = config["FASTQ_DIR"]
include: "functions/BacWork_function.py"
include: "functions/mail.py"
include: "functions/parsing_files.py"

### import rules
include: "rules/reads_preprocessing.smk"
include: "rules/assembly_and_annotation.smk"
include: "rules/sequence_typing.smk"
include: "rules/Bacillus.smk"
include: "rules/Escherichia.smk"
include: "rules/Campylobacter.smk"
include: "rules/Listeria.smk"
include: "rules/Salmonella.smk"
include: "rules/Staphylococcus.smk"
include: "rules/Clostridium.smk"
include: "rules/resume.smk"

makedirs(WORKDIR)

rule all:
	input:
		get_all_input(config)

onsuccess:
	send_mail_onsuccess(config)
	clean_tmp_files()

onerror:
	send_mail_onerror(config)