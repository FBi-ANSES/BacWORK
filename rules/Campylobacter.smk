rule Campylobacter_abricate_host_markers:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/Campylobacter/abricate/host_markers.tsv"
	params:
		mincov = config["abricate"]["mincov"],
		minid = config["abricate"]["minid"],
	conda:
		"../envs/abricate.yaml"
	threads:
		config["abricate"]["threads"]
	shell:"""
		mkdir -p `dirname {output}` 
		abricate --threads {threads} --datadir db/abricate_db --db campy --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output}
		"""		