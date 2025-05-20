configfile: "config/config_plascopedb.json"


rule all:
	input:
		f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb_uniq.fna.nsq",
		f"{config['database_directory']}/{config['genus']}_{config['species']}/seqid_to_taxid.map",
		f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_pls_chr.fna",
		f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plascopedb.1.cf"


rule extract_plasmid_sequences:
	input:
		PLSdb_fasta = config['PLSdb_fasta']
	params:
		genus = config["genus"],
		species = config["species"],
		outdir = config["database_directory"]
	output:
		fasta = temp(f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb.fna")
	shell:
		"""
		mkdir -p {params.outdir}/{params.genus}_{params.species} ;
		python scripts/grepSeq.py -p {input.PLSdb_fasta} -i \"{params.genus}\" -o {output.fasta} ;
		sed -i 's/ .*$//g' {output.fasta}
		"""


rule clustering_plasmid_sequences:
	input:
		fasta = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb.fna"
	params:
		ram = config["ram"]
	output:
		f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb_uniq.fna"
	conda:
		"envs/cdhit.yaml"
	threads:
		config["threads"]
	log:
		"logs/make_plascopedb.log"
	shell:
		"scripts/cd-hit-est -i {input.fasta} -c 1 -M {params.ram} -T {threads} -o {output}  2>> {log}"


rule make_blastdb_nucl:
	input:
		"{fasta}"
	output:
		"{fasta}.nsq"
	conda:
		"envs/blast.yaml"	
	shell:
		"makeblastdb -dbtype nucl -in {input}"


rule download_genomes:
	params:
		genus = config["genus"],
		species = config["species"],
		outdir = config["database_directory"]
	output:
		chr = temp(f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_chr.fna")
	conda:
		"envs/ncbi-genome-download.yaml"
	shell:
		"ncbi-genome-download --genera '{params.genus}' --formats fasta --assembly-levels complete,chromosome bacteria ;\
		gunzip refseq/bacteria/*/*gz ; \
		cat refseq/bacteria/*/*fna > {output.chr} ; \
		rm -r refseq ; \
		sed -i 's/ .*$//g' {output.chr}"


rule remove_plasmids_from_chromosomes:
	input:
		pls = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb_uniq.fna",
		chr = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_chr.fna"
	output:
		temp(f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_chr_uniq.fna")
	shell:
		"python scripts/cleanChr.py -p {input.pls} -c {input.chr} -o {output}"
		
		
rule make_map_file:
	input:
		pls = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plsdb_uniq.fna",
		chr = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_chr_uniq.fna"
	output:
		temp(f"{config['database_directory']}/{config['genus']}_{config['species']}/seqid_to_taxid.map")
	run:
		shell("grep '>' {input.pls} | sed 's/>//g' | sed 's/$/\t3/g' > {output}")
		shell("grep '>' {input.chr} | sed 's/>//g' | sed 's/$/\t2/g' >> {output}")
		
		
rule merge_chr_pls:
	input :
		chr = "{prefix}_chr_uniq.fna",
		pls = "{prefix}_plsdb_uniq.fna"
	output:
		temp("{prefix}_pls_chr.fna")
	shell:
		"cat {input.chr} {input.pls} > {output}"


rule centrifuge_build:
	input:
		fasta = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_pls_chr.fna",
		map = f"{config['database_directory']}/{config['genus']}_{config['species']}/seqid_to_taxid.map"
	params:
		taxonomy_tree = config["taxonomy-tree"],
		name_table = config["name-table"],
		output_prefix = f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plascopedb"
	output:
		f"{config['database_directory']}/{config['genus']}_{config['species']}/all_{config['genus']}_{config['species']}_plascopedb.1.cf"
	threads:
		config["threads"]
	conda:
		"envs/plascope.yaml"
	shell:
		"centrifuge-build -p {threads} --conversion-table {input.map} --taxonomy-tree {params.taxonomy_tree} --name-table {params.name_table} {input.fasta} {params.output_prefix}"
		
		
