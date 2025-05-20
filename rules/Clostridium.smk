rule Clostridium_abricate_virulence_markers:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/Clostridium/abricate/virulence_markers.tsv"
	params:
		mincov = config["abricate"]["mincov"],
		minid = config["abricate"]["minid"],
	conda:
		"../envs/abricate.yaml"
	threads:
		config["abricate"]["threads"]
	shell:"""
		mkdir -p `dirname {output}` 
		abricate --threads {threads} --datadir db/abricate_db --db VCP --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output}
		"""		