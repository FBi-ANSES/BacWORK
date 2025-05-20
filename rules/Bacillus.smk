rule Bacillus_btyper3:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		result = f"{WORKDIR}/{{sample}}/Bacillus/btyper3/{{sample}}_final_results.txt"
	params:
		outdir = f"{WORKDIR}/{{sample}}/Bacillus",
		btyper3_outdir = f"{WORKDIR}/{{sample}}/Bacillus/btyper3_final_results",
		final_dir = f"{WORKDIR}/{{sample}}/Bacillus/btyper3"
	conda:
		"../envs/btyper3.yaml"
	shell:"""
		rm -r {params.final_dir}
		btyper3  --ani_species True --ani_subspecies True --ani_typestrains True --bt True --mlst True --panC True --download_mlst_latest True -i {input.assembly} -o {params.outdir}
		mkdir {params.final_dir}
		mv {params.btyper3_outdir}/* {params.final_dir}
		rm -r {params.btyper3_outdir}
		"""

rule Bacillus_Bt_detect:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		result = f"{WORKDIR}/{{sample}}/Bacillus/Bt_detect/{{sample}}.tsv",
		blast = temp(f"{{sample}}_blast.tsv")
	conda:
		"../envs/virulencefinder.yaml"
	threads:
		config["blast"]["threads"]		
	shell:"""
		python scripts/Bt_detect.py -t db/table_bt.txt -db db/abricate_db/ID_Bt/sequences -T {threads} -i {input.assembly} > {output.result}
		"""
		
rule Bacillus_abricate_specific_virulence:
	input:
		assembly = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa"
	output:
		f"{WORKDIR}/{{sample}}/Bacillus/abricate/VBC.tsv"
	params:
		mincov = config["abricate"]["mincov"],
		minid = config["abricate"]["minid"],
	conda:
		"../envs/abricate.yaml"
	threads:
		config["abricate"]["threads"]
	shell:"""
		mkdir -p `dirname {output}` 
		abricate --threads {threads} --datadir db/abricate_db --db VBC --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output}
		"""				
