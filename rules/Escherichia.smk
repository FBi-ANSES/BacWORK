rule Escherichia_clermontyping:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	params:
		dir = f"{WORKDIR}/{{sample}}/Escherichia/clermontyping",
		out = f"{WORKDIR}/{{sample}}/Escherichia/clermontyping/analysis_*/analysis_*phylogroups.txt"
	output:
		f"{WORKDIR}/{{sample}}/Escherichia/clermontyping/{{sample}}_phylogroups.txt"
	conda:
		"../envs/clermontyping.yaml"
	shell:
		"""
		rm -r {params.dir} ;
		mkdir -p {params.dir} ;
		cp -r scripts/ClermonTyping {params.dir}/. ;
		touch {output} ; 
		cd {params.dir} ; 
		ClermonTyping/clermonTyping.sh --fasta {input}  ; 
		cd ../../../.. ; 
		cp {params.out} {output}
		"""
				
rule Escherichia_stecfinder:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		result = f"{WORKDIR}/{{sample}}/Escherichia/stecfinder/{{sample}}_stecfinder.txt"
	conda:
		"../envs/stecfinder.yaml"
	threads:
		config["stecfinder"]["threads"]
	shell:
		"""
		stecfinder -i {input} --hits --output {output} -t {threads}
		"""				