
if config["umi_tools"]["extract"]["five_prime"]:
	umi_extract_extra = config["umi_tools"]["extract"]["extra"]+""
else:
	umi_extract_extra = config["umi_tools"]["extract"]["extra"]+" --3prime"


if config["umi_tools"]["trimm"]["trimm"]:
	rule umi_extract:
		input:
			fastq=IN+"/{sample}"+config["fastq"]
		output:
			fastq=OUT+"/umi/{sample}"+config["fastq"]
		log:OUT+"/umi/{sample}.log"
		conda:
			"../envs/umi_and_trim.yaml"
		params:
			trimfile=OUT+"/umi/{sample}_fastx"+config["fastq"],
			trimnumber=config["umi_tools"]["trimm"]["trim_supp"]+1,
			extra= umi_extract_extra,
			bc_pattern=config["umi_tools"]["extract"]["bc-pattern"],
			extract_method=config["umi_tools"]["extract"]["extract-method"]
		shell:
			"fastx_trimmer -f {params.trimnumber} -i {input.fastq} -o {params.trimfile} && "
			"umi_tools extract --bc-pattern={params.bc_pattern} --extract-method={params.extract_method} --log={log} --stdout={output.fastq} --stdin={params.trimfile} {params.extra} && "
			"rm {params.trimfile}"
	
else:
	rule umi_extract:
		input:
			fastq=IN+"/{sample}"+config["fastq"]
		output:
			fastq=OUT+"/umi/{sample}"+config["fastq"]
		log:OUT+"/umi/{sample}.log"
		conda:
			"../envs/umi_tools.yaml"
		params:
			extra= umi_extract_extra,
			bc_pattern=config["umi_tools"]["extract"]["bc-pattern"],
			extract_method=config["umi_tools"]["extract"]["extract-method"]
		shell:
			"umi_tools extract --bc-pattern={params.bc_pattern} --extract-method={params.extract_method} --log={log} --stdout={output.fastq} --stdin={input.fastq} {params.extra}"

rule umi_dedup:
	input:
		bam=OUT+"/{genome}/mapping/bam/sorted/{prefix}.sorted.bam",
		bai=OUT+"/{genome}/mapping/bam/sorted/{prefix}.sorted.bai"
	output:
		bam=OUT+"/{genome}/mapping/bam/umidedup/{prefix}.umidedup.bam"
	conda:
		"../envs/umi_tools.yaml"
	params:
		OUT+"/{genome}/mapping/bam/umidedup/{prefix}"
	shell:
		"umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output.bam}"
