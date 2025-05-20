rule Staphylococcus_spatyper:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		spatyper = f"{WORKDIR}/{{sample}}/Staphylococcus/spatyper/{{sample}}_spatyper.txt"
	params:
		sparepeat = config["spaTyper"]["spa_repeat_file"],
		repeatorder = config["spaTyper"]["repeat_order_file"],
		folder = f"{WORKDIR}/{{sample}}/Staphylococcus/spatyper"
	conda:
		"../envs/spatyper.yaml"
	shell:
		"""
		spaTyper -r {params.sparepeat} -o {params.repeatorder} -f {input.assembly} -d {params.folder} --output {output.spatyper}
		"""
		
rule Staphylococcus_staphopia_sccmec:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/Staphylococcus/staphopia-sccmec/{{sample}}_staphopia-sccmec.txt"
	conda:
		"../envs/staphopia-sccmec.yaml"
	shell:"""
		 staphopia-sccmec --assembly {input.assembly} > {output}
		"""		
		
rule Staphylococcus_naura:
	input:
		gbk = rules.assembly_and_annotation_bakta.output.gbk
	output:
		output = f"{WORKDIR}/{{sample}}/Staphylococcus/naura/matrix.tsv"
	params:
		others = config["naura"]["others"]
	conda: 
		"../envs/naura.yaml"
	threads:
		config["naura"]["threads"]
	shell:
		"""
		cp {input.gbk} `dirname {output.output}`/bakta.gbk
		cp -r db/NAuRA_staph_prot `dirname {output.output}`/.
		cp db/NAuRA_staph_prot_list.txt `dirname {output.output}`/.
		cd `dirname {output.output}`
		sed "s|^|$(pwd)/|" NAuRA_staph_prot_list.txt
		NAuRA -i . \
		-q NAuRA_staph_prot_list.txt \
		-T {threads} \
		{params.others}
		"""
		