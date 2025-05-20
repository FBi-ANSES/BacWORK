rule sequence_typing_mlst:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		mlst_report = f"{WORKDIR}/{{sample}}/sequence_typing/mlst/{{sample}}_mlst.txt"
	params:
		novel_allele_file = f"{WORKDIR}/{{sample}}/mlst/{{sample}}_novel.fa",
	conda:
		"../envs/mlst.yaml"
	threads:
		config["mlst"]["threads"]
	shell:
		"""
		mlst --threads {threads} --nopath --novel {params.novel_allele_file} \
		--label {wildcards.sample} {input} > {output}
		"""

rule sequence_typing_chewbbacaSchema:
	input:
		fasta = get_assembly_list(config)
	output:
		schema = f"{WORKDIR}/cgMLST/.end"
	params:
		config_json = config['samples'],
		schema = config["cgmlst_database"]
	conda:
		"../envs/chewbbaca.yaml"
	threads:
		config["chewbbaca"]["threads"]
	script:
		"../scripts/chewbbacaSchema.py"

rule sequence_typing_chewbbaca:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		schema = f"{WORKDIR}/cgMLST/.end"
	output:
		profile = f"{WORKDIR}/{{sample}}/sequence_typing/chewbbaca/alleles.tsv"
	params: 
		genus = lambda wildcards:  get_genus(wildcards.sample),
	conda:
		"../envs/chewbbaca.yaml"
	script:
		"../scripts/chewbbaca.py"
        
rule sequence_typing_chewbbacaHash:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		schema = f"{WORKDIR}/cgMLST/.end"
	output:
		profile = f"{WORKDIR}/{{sample}}/sequence_typing/chewbbaca/alleles_hash.tsv"
	conda:
		"../envs/chewbbaca.yaml"		
	params: 
		genus = lambda wildcards:  get_genus(wildcards.sample),
	script:
		"../scripts/chewbbacaHash.py"
        

# rule sequence_typing_chewbbaca:
	# input:
		# genome_list = f"{WORKDIR}/cgMLST/genome_list.txt",
	# output:
		# allele = f"{WORKDIR}/cgMLST/alleles.tsv",
		# contigsInfo = f"{WORKDIR}/cgMLST/contigsInfo.tsv",
		# statistics = f"{WORKDIR}/cgMLST/statistics.tsv"
	# threads: config["chewbbaca"]["threads"]
	# conda:"../envs/chewbbaca.yaml"
	# params:
		# schema = config["cgmlst_database"]+config["chewbbaca"][GENUS]
	# shell:
		# """ cp {params.schema} `dirname {output.allele}`/schema/
			# chewBBACA.py AlleleCall -i {input.genome_list} -g `dirname {output.allele}`/schema/ -o `dirname {output.allele}` --cpu {threads} --gl `dirname {output.allele}`/schema/listGenes.txt && 
			# cp `dirname {output.allele}`/results_*/results_alleles.tsv `dirname {output.allele}`/alleles.tsv && 
			# cp `dirname {output.allele}`/results_*/results_contigsInfo.tsv `dirname {output.allele}`/contigsInfo.tsv && 
			# cp `dirname {output.allele}`/results_*/results_statistics.tsv `dirname {output.allele}`/statistics.tsv && 
			# rm -r `dirname {output.allele}`/results_*"""


# rule sequence_typing_chewbbaca_quality:
	# input:
	   # allele = "{run}/cgMLST/alleles.tsv"
	# output:
		# quality = "{run}/cgMLST/quality.end"
	# conda: "../envs/chewbbaca.yaml"
	# shell:
		# """chewBBACA.py TestGenomeQuality -i {input.allele} -n 12 -t 200 -s 5 -o `dirname {output.quality}` &&
			# touch {output.quality}"""
		
		
# rule sequence_typing_chewbbaca_sample:
	# input:
		# fasta = "{run}/{smp}/assembly/shovill/{smp}_contigs_filter.fa",
		# allele = "{run}/cgMLST/alleles.tsv",
		# quality = "{run}/cgMLST/quality.end"
	# output:
		# allele = "{run}/{smp}/cgMLST/chewbbaca_allele/alleles.tsv"
	# threads:1
	# script:
		# "../scripts/chewbbaca_alleleSample.py"

# rule sequence_typing_chewbbaca_hash:
	# input:
		# fasta = "{run}/{smp}/assembly/shovill/{smp}_contigs_filter.fa",
		# contigsInfo = "{run}/cgMLST/contigsInfo.tsv",
		# quality = "{run}/cgMLST/quality.end"
	# output:
		# allele = "{run}/{smp}/cgMLST/chewbbaca_allele/alleles_hash.tsv"
	# threads:1
	# script:
		# "../scripts/chewbbaca_hash.py"