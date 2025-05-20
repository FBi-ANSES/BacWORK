rule Listeria_CC:
	input:
		mlst_report = f"{WORKDIR}/{{sample}}/sequence_typing/mlst/{{sample}}_mlst.txt",
		CC = "db/Listeria_CC.tsv"
	output:
		CC_report = f"{WORKDIR}/{{sample}}/Listeria/sequence_typing/mlst/{{sample}}_CC.txt"
	script:
		"../scripts/Listeria_CC.py"

rule Listeria_serotyping:
	input:
		contigs = f"{WORKDIR}/{{sample}}/assembly_and_annotation/filter/{{sample}}.fa",
		amplicon = "db/Listeria_amplicon.fa"
	output:
		blast = f"{WORKDIR}/{{sample}}/Listeria/serotyping/{{sample}}_blastn.txt"
	conda: 
		"../envs/seqsero2.yaml"
	shell:
		"""
		makeblastdb -dbtype nucl -out `dirname {output.blast}`/contigs.fa -in {input.contigs}
		blastn -db `dirname {output.blast}`/contigs.fa -query {input.amplicon} \
		-out {output.blast} -outfmt "6 qseqid sseqid evalue sstart send pident"
		"""


rule Listeria_serotyping_report:
	input:
		blast = f"{WORKDIR}/{{sample}}/Listeria/serotyping/{{sample}}_blastn.txt"
	output:
		serotype = f"{WORKDIR}/{{sample}}/Listeria/serotyping/{{sample}}_serotyping.txt"
	script:
		"../scripts/Listeria_serotyping_report.py"
