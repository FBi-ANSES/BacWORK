rule assembly_and_annotation_shovill :
	input:
		fastq_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		fastq_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz"
	output :
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/shovill/contigs.fa",
	params:
		tmp = f"{WORKDIR}/{{sample}}/assembly_and_annotation/shovill/tmp/",
		minlen = config["shovill"]["minlen"],
		ram = config["shovill"]["ram"],
		depth = config["shovill"]["depth"],
	threads: 
		config["shovill"]["threads"]
	conda: 
		"../envs/shovill.yaml"
	shell:
		"""
		shovill --R1 {input.fastq_R1} --R2 {input.fastq_R2} --outdir `dirname {output.contigs}` --force \
		--minlen {params.minlen} \
		--cpus {threads} \
		--tmpdir {params.tmp} \
		--depth {params.depth} \
		--opts "-m {params.ram} --tmp-dir {params.tmp}/spades" \
		--ram {params.ram}
		"""

rule assembly_and_annotation_mauve_reorder:
	input:
		shovill_assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/shovill/contigs.fa",
		ref_fasta = rules.reads_preprocessing_download_reference.output.ref_fasta
	output:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/mauvecm/contigs.fa.fas",
	params:
		res_dir = f"{WORKDIR}/{{sample}}/assembly_and_annotation/mauvecm"
	conda:
		"../envs/mauvecm.yaml"
	threads:
		config["mauvecm"]["threads"]
	shell:
		"""
		perl scripts/NGS_MauveReorder_on_ref.pl -fa_ref {input.ref_fasta} -draft {input.shovill_assembly} -p {threads} -res_dir {params.res_dir}
		"""

rule assembly_and_annotation_resfinder:
	input:
		R1 = rules.reads_preprocessing_fastp.output.fastq_R1, 
		R2 = rules.reads_preprocessing_fastp.output.fastq_R2, 
		kraken2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/kraken2/kraken_report.csv"
	output:
		results = f"{WORKDIR}/{{sample}}/assembly_and_annotation/resfinder/ResFinder_results.txt",
		version = f"{WORKDIR}/{{sample}}/assembly_and_annotation/resfinder/ResFinder_version.txt"
	params:
		resfinder_db = config["resfinder_database"],
		pointfinder_db = config["pointfinder_database"],
		disinfinder_db = config["disinfinder_database"],
		species = lambda wildcards, input:  get_species(wildcards.sample,input.kraken2),
		output_dir = f"{WORKDIR}/{{sample}}/assembly_and_annotation/resfinder"
	conda:
		"../envs/resfinder.yaml"
	shell:"""
		python -m resfinder -ifq {input.R1} {input.R2} -o {params.output_dir} \
		-s {params.species} --ignore_missing_species -acq -d -c \
		-db_res {params.resfinder_db} \
		-db_disinf {params.disinfinder_db} \
		-db_point {params.pointfinder_db}
		python -m resfinder -v > {output.version}
	"""

rule assembly_and_annotation_virulencefinder:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/assembly_and_annotation/virulencefinder/results_tab.tsv"
	params:
		virulencefinder_db = config["virulencefinder_database"]
	conda:
		"../envs/virulencefinder.yaml"
	shell:"""
		python scripts/virulencefinder.py -i {input.assembly} -o `dirname {output}` -x -p {params.virulencefinder_db}
		"""

rule assembly_and_annotation_filter :
	"""
	Filtre des contigs sur: mincov
	"""
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/mauvecm/contigs.fa.fas",
		plascope_dir = directory(f"{WORKDIR}/{{sample}}/assembly_and_annotation/plascope/{{sample}}_PlaScope/PlaScope_predictions")
	output :
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		report = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/report.tsv"
	params:
		mincov = config["assembly_and_annotation_filter"]["mincov"],
		sample = "{sample}"
	script:
		'../scripts/assembly_and_annotation_filter.py'

rule assembly_and_annotation_quastReads :
	"""
	Qualité de l'assemblage avec quast
	Calcul du nombre de contigs, taille cumulée de l'assemblage et N50
	"""
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		fastq_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		fastq_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
		ref_fasta = rules.reads_preprocessing_download_reference.output.ref_fasta
	output :
		quast = f"{WORKDIR}/{{sample}}/assembly_and_annotation/quast/report.tsv",
	threads: config['quast']['threads']	
	conda:
		"../envs/quast.yaml"
	shell:
		"""
		quast -o `dirname {output.quast}` --pe1 {input.fastq_R1} --pe2 {input.fastq_R2} \
		-t {threads} -r {input.ref_fasta} {input.contigs}
		"""

rule assembly_and_annotation_bakta:
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		gbk=f"{WORKDIR}/{{sample}}/assembly_and_annotation/bakta/{{sample}}.gbff",
		gff=f"{WORKDIR}/{{sample}}/assembly_and_annotation/bakta/{{sample}}.gff3",
	threads: config['bakta']['threads']
	params: 
		db=config['bakta_database']
	conda:
		"../envs/bakta.yaml"
	shell:
		"""
		bakta {input.contigs} --db {params.db} --output `dirname {output.gbk}` --threads {threads} --force
		"""

rule assembly_and_annotation_busco: 
	input: 
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		busco = f"{WORKDIR}/{{sample}}/assembly_and_annotation/busco/short_summary.json",
		tmp_contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/busco/{{sample}}.fa",
	threads: config['busco']['threads']
	params:
		mode = config['busco']['mode'],
		lineage = config['busco']['lineage'],
		db = config['busco_database'],
		other_options = config['busco']['other_options'],
		workdir = f"{WORKDIR}",
		outdir = f"{{sample}}/assembly_and_annotation/busco/"
	conda:
		"../envs/busco.yaml"
	shell:
		"""
		echo "Checking variables:"
		echo "smp: {wildcards.sample}"
		cp {input.contigs} {output.tmp_contigs}
		cd `dirname {output.busco}`
		busco -i {output.tmp_contigs} -o analyse_busco -m {params.mode} -f --cpu {threads} --download_path {params.db} {params.other_options} {params.lineage}
		cp analyse_busco/short_summary.generic*.json short_summary.json
		"""

rule assembly_and_annotation_abricate_vfdb:
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa" 
	output:
		f"{WORKDIR}/{{sample}}/assembly_and_annotation/abricate/vfdb.tsv"
	params:
		mincov = config["abricate"]["mincov"],
		minid = config["abricate"]["minid"],
	conda:
		"../envs/abricate.yaml"
	threads:
		config["abricate"]["threads"]
	shell:"""
		mkdir -p `dirname {output}` 
		abricate --threads {threads} --db vfdb --minid {params.minid} --mincov {params.mincov} {input.contigs} > {output}
		"""

rule assembly_and_annotation_vibrant :
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/assembly_and_annotation/vibrant/VIBRANT_phages.fa"
	params:
		databases_dir = config["vibrant_database"],
		file_dir = config["vibrant_files"],
		out_dir = f"{WORKDIR}/{{sample}}/assembly_and_annotation/vibrant/"
	conda:
		"../envs/vibrant.yaml"
	threads:
		config["vibrant"]["threads"]
	shell:"""
		rm -r `dirname {output}`
		VIBRANT_run.py -i {input.contigs} -f nucl -folder {params.out_dir} \
		-d {params.databases_dir} -m {params.file_dir} -t {threads}
		cp {params.out_dir}VIBRANT_*/VIBRANT_phages_*/*.phages_combined.fna {output}
		"""

rule assembly_and_annotation_plasmidfinder:
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/assembly_and_annotation/plasmidfinder/data.json"
	params:
		plasmidfinder_db = config["plasmidfinder_database"]	
	conda:
		"../envs/plasmidfinder.yaml"
	shell:"""
		mkdir -p `dirname {output}` 
		plasmidfinder.py -i {input.contigs} -o `dirname {output}` -p {params.plasmidfinder_db}
		"""

rule assembly_and_annotation_plascope :
        input:
                assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/mauvecm/contigs.fa.fas",
                kraken2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/kraken2/kraken_report.csv"
        params:
                database_dir = lambda wildcards, input:  get_database_dir(wildcards.sample,config['plascope_database'],input.kraken2),
                database_name = lambda wildcards, input:  get_database_name(wildcards.sample,config['plascope_database'],input.kraken2),
                contigs_rename = f"{WORKDIR}/{{sample}}/assembly_and_annotation/plascope/contigs_rename.fa",
                output_dir = f"{WORKDIR}/{{sample}}/assembly_and_annotation/plascope"
        threads :
                config['plascope']['threads']
        output:
                directory(f"{WORKDIR}/{{sample}}/assembly_and_annotation/plascope/{{sample}}_PlaScope/PlaScope_predictions"),
        conda:
                "../envs/plascope.yaml"
        shell:"""
                rm -r {params.output_dir} ; mkdir {params.output_dir} ;
                cat {input.assembly} | perl -p -e 's/^.* spades=/>/' | sed 's/>.*origname=/>/g' | sed 's/ .*//g' > {params.contigs_rename} ;
                plaScope.sh --fasta {params.contigs_rename} -o {params.output_dir} --db_dir {params.database_dir} --db_name {params.database_name} --sample {wildcards.sample}
                """
