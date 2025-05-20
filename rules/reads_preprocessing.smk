rule reads_preprocessing_fastqc_raw:
	input:
		R1 = f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz", 
		R2 = f"{FASTQ_DIR}/{{sample}}_R2.fastq.gz"
	output:
		html_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastqc/raw/{{sample}}_R1_fastqc.html",
		html_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastqc/raw/{{sample}}_R2_fastqc.html"
	conda:
		"../envs/fastqc.yaml"
	threads:
		config["fastqc"]["threads"]
	shell:
		"""
		fastqc -o `dirname {output.html_R1}` -t {threads} {input.R1} {input.R2} 
		"""


rule reads_preprocessing_fastp:
	input:
		R1 = f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz", 
		R2 = f"{FASTQ_DIR}/{{sample}}_R2.fastq.gz"
	output:
		fastq_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		fastq_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
		json = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}.json",
		html = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}.html"
	params:
		tmp_dir = f"{WORKDIR}/{{sample}}",
		options = config["fastp"]["options"]
	conda:
		"../envs/fastp.yaml"
	threads:
		config["fastp"]["threads"]
	shell:
		"""
		fastp -i {input.R1} -o {output.fastq_R1} -I {input.R2} -O {output.fastq_R2} -j {output.json} -h {output.html} -w {threads} {params.options}
		"""


rule reads_preprocessing_fastqc_clean:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
	output:
		html_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastqc/clean/{{sample}}_clean_1_fastqc.html",
		html_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastqc/clean/{{sample}}_clean_2_fastqc.html"
	conda:
		"../envs/fastqc.yaml"
	threads:
		config["fastqc"]["threads"]
	shell:
		"""
		fastqc -o `dirname {output.html_R1}` -t {threads} {input.R1} {input.R2} 
		"""


rule reads_preprocessing_confindr:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz", 
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz" 
	output:
		report = f"{WORKDIR}/{{sample}}/reads_preprocessing/confindr/confindr_report.csv",
	conda:
		"../envs/confindr.yaml"
	params:
		db = config["confindr_database"]
	threads:
		config["confindr"]["threads"]
	shell:
		"""
		confindr.py -i `dirname {input.R1}` -o `dirname {output.report}` -fid _clean_1 -rid _clean_2 -d {params.db} --rmlst --threads {threads}
		"""


rule reads_preprocessing_kraken2:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
	output:
		report = f"{WORKDIR}/{{sample}}/reads_preprocessing/kraken2/kraken_report.csv",
		kreport = f"{WORKDIR}/{{sample}}/reads_preprocessing/kraken2/out.krepport"
	conda:
		"../envs/kraken2.yaml"
	params:
		db = config["kraken2_database"],
		confidence = config["kraken2"]["confidence"]
	threads:
		config["kraken2"]["threads"]
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --confidence {params.confidence} --paired --report {output.report} --output {output.kreport} {input.R1} {input.R2}
		"""

rule reads_preprocessing_contamination:
	input:
		kraken2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/kraken2/kraken_report.csv",
		confindr = f"{WORKDIR}/{{sample}}/reads_preprocessing/confindr/confindr_report.csv"
	output:
		report = f"{WORKDIR}/{{sample}}/reads_preprocessing/contamination/report.tsv"
	params:
		genus = lambda wildcards:  get_genus(wildcards.sample),
		species = lambda wildcards, input:  get_species(wildcards.sample,input.kraken2),
		confindr_SNVs_threshold = config["confindr"]["SNVs_threshold"]
	shell:
		"""
		python scripts/reads_preprocessing_contamination.py -c {input.confindr} \
		-k {input.kraken2} -g {params.genus} -s {params.species} -t {params.confindr_SNVs_threshold} \
		-o {output.report}
		"""


rule reads_preprocessing_mash_screen:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
	params:
		pvalue = config["mash"]["pvalue"],
		exclude = config["mash"]["exclude"],
		sketch = config["complete_genome_sketch_database"],
		genus = lambda wildcards:  get_genus(wildcards.sample)
	output:
		csv = f"{WORKDIR}/{{sample}}/reads_preprocessing/mash_screen/closer_references.csv",
	conda:
		"../envs/mash.yaml"
	threads:
		config["mash"]["threads"]
	shell:
		"""
		mash screen -v {params.pvalue} -p {threads} {params.sketch} {input.R1} {input.R2} |\
		grep -v {params.exclude} |\
		grep {params.genus} |\
		sort -r -k 1 > {output.csv}
		"""


rule reads_preprocessing_download_reference:
	input:
		mash_screen_output = rules.reads_preprocessing_mash_screen.output.csv
	params:
		accession = lambda wildcards, input: grep_accession_from_mash_screen_output(input.mash_screen_output)
	output:
		ref_fasta = f"{WORKDIR}/{{sample}}/reads_preprocessing/mash_screen/reference.fa",
		ref_gff = f"{WORKDIR}/{{sample}}/reads_preprocessing/mash_screen/reference.gff"
	conda:
		"../envs/ncbi-acc-download.yaml"
	shell:
		"""
		cd `dirname {input}` ;
		ncbi-acc-download --format fasta {params.accession} ;
		ncbi-acc-download --format gff3 {params.accession} ;
		mv {params.accession}.fa reference.fa ;
		mv {params.accession}.gff reference.gff ;
		"""


rule reads_preprocessing_bbmap:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
		ref_fasta = rules.reads_preprocessing_download_reference.output.ref_fasta
	params:
		k = config["bbmap"]["k"],
		covbinsize = config["bbmap"]["covbinsize"],
		other = config["bbmap"]["other_options"],
		ram = config["bbmap"]["ram"]
	output:
		covstats = f"{WORKDIR}/{{sample}}/reads_preprocessing/bbmap/covstats.txt"
	conda:
		"../envs/bbmap.yaml"
	threads:
		config["bbmap"]["threads"]
	shell:
		"""
		bbmap.sh in={input.R1} in2={input.R2} ref={input.ref_fasta} \
		covstats={output.covstats} covbinsize={params.covbinsize} \
		k={params.k} threads={threads} \
		-Xmx{params.ram} {params.other}
		"""	


rule reads_preprocessing_speciesfinder:
	input:
		R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz",
	params:
		db = config["speciesfinder_database"],
		outdir = f"{WORKDIR}/{{sample}}/reads_preprocessing/speciesFinder/",
	output:
		json = f"{WORKDIR}/{{sample}}/reads_preprocessing/speciesFinder/data.txt" 
	conda:
		"../envs/cge.yaml"
	shell:
		"""
		python scripts/speciesfinder.py -i {input.R1} {input.R2} \
		-p {params.db} -o {params.outdir} -tmp {params.outdir}
		"""

