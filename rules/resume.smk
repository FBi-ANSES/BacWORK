rule resume_make_summary_json:
	input:
		fastp_output = rules.reads_preprocessing_fastp.output.json,
		clean_R1 = rules.reads_preprocessing_fastp.output.fastq_R1,
		bbmap_covstat = rules.reads_preprocessing_bbmap.output.covstats,
		quast_report = rules.assembly_and_annotation_quastReads.output.quast,
		kraken_report = rules.reads_preprocessing_kraken2.output.report,
		mlst_report = rules.sequence_typing_mlst.output.mlst_report,
		assembly = rules.assembly_and_annotation_filter.output.contigs,
		contamination_report = rules.reads_preprocessing_contamination.output.report,
	output:
		json = f"{WORKDIR}/{{sample}}/resume/summary.json"
	params:
		sample_id = f"{{sample}}"
	run:
		make_resume_json(params.sample_id,output.json,input.fastp_output,\
		input.clean_R1,input.bbmap_covstat,input.quast_report,input.kraken_report,\
		input.mlst_report,input.assembly,input.contamination_report)

rule resume_make_tools_list:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		resfinder_version = rules.assembly_and_annotation_resfinder.output.version,
	output:
		outfile = f"{WORKDIR}/{{sample}}/resume/tools_list.tsv"
	run:
		make_tools_list(output.outfile,input.resfinder_version)
		
rule resume_make_db_version_file:
	output:
		outfile = f"{WORKDIR}/{{sample}}/resume/db_version.tsv"
	run:
		make_db_version_file(output.outfile)
		
rule resume_copy_quality_json:
	input:
		"config/config_quality.json"
	output:
		f"{WORKDIR}/{{sample}}/resume/quality_threshold.json"
	shell:
		"""
		 cp {input} {output}
		"""		
		
rule resume_copy_config_tools_json:
	input:
		"config/config_tools.json"
	output:
		f"{WORKDIR}/{{sample}}/resume/config_tools.json"
	shell:
		"""
		 cp {input} {output}
		"""
		
rule resume_copy_assembly:
	input:
		assembly = rules.assembly_and_annotation_filter.output.contigs
	output:
		f"{WORKDIR}/{{sample}}/resume/{{sample}}_assembly.fa"
	shell:
		"""
		 cp {input} {output}
		"""		
		
rule resume_suppl_Salmonella:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		seqsero_tsv = rules.Salmonella_seqsero2.output.tsv,
		sistr_tsv = rules.Salmonella_sistr.output.sister_output
	output:
		flag = f"{WORKDIR}/{{sample}}/resume/resume_Salmonella"
	run:
		add_Salmo_informations(input.resume_json, input.seqsero_tsv, input.sistr_tsv)
		os.system(f"touch {output.flag}")
		
rule resume_suppl_Bacillus:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		bt_detect = rules.Bacillus_Bt_detect.output.result,
		btyper = rules.Bacillus_btyper3.output.result
	output:
		flag = f"{WORKDIR}/{{sample}}/resume/resume_Bacillus"
	run:
		add_Bacillus_informations(input.resume_json, input.bt_detect, input.btyper)
		os.system(f"touch {output.flag}")

rule resume_suppl_Staphylococcus:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		spatyper = rules.Staphylococcus_spatyper.output.spatyper,
	output:
		flag = f"{WORKDIR}/{{sample}}/resume/resume_Staphylococcus"
	run:
		add_Staph_informations(input.resume_json, input.spatyper)
		os.system(f"touch {output.flag}")
		
rule resume_suppl_Escherichia:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		stecfinder = rules.Escherichia_stecfinder.output.result,
	output:
		flag = f"{WORKDIR}/{{sample}}/resume/resume_Escherichia"
	run:
		add_Escherichia_informations(input.resume_json, input.stecfinder)
		os.system(f"touch {output.flag}")
		
rule resume_suppl_Listeria:
	input:
		resume_json = rules.resume_make_summary_json.output.json,
		serotyping_report = rules.Listeria_serotyping_report.output.serotype,
		CC_report = rules.Listeria_CC.output.CC_report,
	output:
		flag = f"{WORKDIR}/{{sample}}/resume/resume_Listeria"
	run:
		add_Listeria_informations(input.resume_json, input.serotyping_report, input.CC_report)
		os.system(f"touch {output.flag}")

rule resume_multiqc:
	input:
		multiqc_input_files(config)
	output:
		html = f"{WORKDIR}/Bacwork_reports/multiqc_report.html"
	params:
		outdir = f"{WORKDIR}/Bacwork_reports",
		fastp = f"{WORKDIR}/*/reads_preprocessing/fastp",
		fastqc = f"{WORKDIR}/*/reads_preprocessing/fastqc",
		bakta = f"{WORKDIR}/*/assembly_and_annotation/bakta",
		quast = f"{WORKDIR}/*/assembly_and_annotation/quast"
	conda:
		"../envs/multiqc.yaml"	
	shell:
		"""
		 multiqc -f --outdir {params.outdir} --clean-up {params.fastp} {params.fastqc} {params.bakta} {params.quast}
		"""
		
rule resume_merge_json:
	input:
		json_list = resume_files_list(config)
	output:
		json = f"{WORKDIR}/Bacwork_reports/resume.json"
	run:
		merge_summary_json(input.json_list, output.json)
		
rule resume_json_to_tab:
	input:
		json = rules.resume_merge_json.output.json
	output:
		tab = f"{WORKDIR}/Bacwork_reports/resume.tsv"
	run:
		json_to_tab(input.json, output.tab)