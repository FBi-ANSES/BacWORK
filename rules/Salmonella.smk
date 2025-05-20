rule Salmonella_sistr:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		sister_output = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_sistr.json",
		cgmlst_output = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_cgmlst.csv",
		new_alleles = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_new_alleles.fa",
		alleles = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_alleles.json"
	params:
		TMP_dir = f"{WORKDIR}/{{sample}}/Salmonella/sistr",
		sister_output = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_sistr"
	conda:
		"../envs/sistr.yaml"
	threads:
		config["sistr"]["threads"]
	shell:
		"""
		sistr -i {input.assembly} {wildcards.sample} -f json -o {params.sister_output} -MM --more-results \
        --run-mash -p {output.cgmlst_output} -n {output.new_alleles} -a {output.alleles} --qc \
		-t {threads} -T {params.TMP_dir}
		"""		


rule Salmonella_seqsero2:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		txt = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2/SeqSero_result.txt",
		tsv = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2/SeqSero_result.tsv"
	params: 
		b = config["seqsero2"]["mapping"],
		m = config["seqsero2"]["modele"]
	conda:
		"../envs/seqsero2.yaml"
	threads:
		config["seqsero2"]["threads"]
	shell:
		"""
		SeqSero2_package.py -i {input.assembly} -t 4 -b {params.b} \
		-p {threads} -m {params.m} -d `dirname {output.txt}`
		"""

rule Salmonella_seqsero2Reads:
	input:
		fastq_R1 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_1.fastq.gz",
		fastq_R2 = f"{WORKDIR}/{{sample}}/reads_preprocessing/fastp/{{sample}}_clean_2.fastq.gz"
	output:
		txt = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2Reads/SeqSero_result.txt",
		tsv = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2Reads/SeqSero_result.tsv"
	conda: 
		"../envs/seqsero2.yaml"
	params:
		b = config["seqsero2"]["mapping"],
		m = config["seqsero2"]["modele"],
		R1 = f"{{sample}}_clean_1.fastq",
		R2 = f"{{sample}}_clean_2.fastq",
		dir = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2Reads/"
	threads:
		config["seqsero2"]["threads"]
	shell:
		"""
		SeqSero2_package.py -m {params.m} -t 2 -b {params.b} \
		-i {input.fastq_R1} {input.fastq_R2} \
		-p {threads} -d {params.dir}
		"""

rule Salmonella_isPCr:
	input:
		contigs =f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		query = "db/Salmonella_typing/query_file"
	output:
		PCR = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_IsPCR.txt"
	conda: 
		"../envs/isPCR.yaml"
	shell:
		"""
		isPcr {input.contigs} {input.query} {output.PCR}
		"""

rule Salmonella_fliA_fliB:
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		query = "db/Salmonella_typing/fliA_fliB_amplicon.fasta"
	output:
		blast = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_fliA_fliB_blastn.txt"
	conda: 
		"../envs/isPCR.yaml"
	shell:
		"""
		makeblastdb -dbtype nucl -out `dirname {output.blast}`/contigs.fa -in {input.contigs}
		blastn -db `dirname {output.blast}`/contigs.fa -query {input.query} -out {output.blast} -outfmt "6 qseqid sseqid evalue sstart send pident"
		"""

rule Salmonella_typhymurium:
	input:
		isPCR_result = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_IsPCR.txt",
		blast_result = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_fliA_fliB_blastn.txt"
	output:
		typhymurium = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_typhimurium_variant.txt"
	conda: 
		"../envs/cge.yaml"	
	script:
		'../scripts/Salmonella_typhimurium.py'

rule Salmonella_serotyping:
	"""
	Génération du rapport json par échantillon, contenant toutes les informations d'intérêt pour produire le rapport final
	Annotation de qualité de l'échantillon
	"""
	input:
		sistr = f"{WORKDIR}/{{sample}}/Salmonella/sistr/{{sample}}_sistr.json",
		seqsero = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2/SeqSero_result.tsv",
		seqseroReads = f"{WORKDIR}/{{sample}}/Salmonella/seqsero2Reads/SeqSero_result.tsv",
		typhymurium = f"{WORKDIR}/{{sample}}/Salmonella/typhivar/{{sample}}_typhimurium_variant.txt"
	output:
		txt = f"{WORKDIR}/{{sample}}/Salmonella/serotyping/{{sample}}.txt",
		serotyping = f"{WORKDIR}/{{sample}}/Salmonella/serotyping/{{sample}}.serotype"
	conda: 
		"../envs/cge.yaml"			
	script:
		'../scripts/Salmonella_serotyping.py'
