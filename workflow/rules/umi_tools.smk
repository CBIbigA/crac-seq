
if config["umi_tools"]["extract"]["five_prime"]:
	umi_extract_extra = config["umi_tools"]["extract"]["extra"]+""
else:
	umi_extract_extra = config["umi_tools"]["extract"]["extra"]+" --3prime"


if config["umi_tools"]["trimm"]["trimm"]:
	rule umi_extract:
		input:
			fastq=IN+"/{sample}"+EXT
		output:
			fastq=OUT+"/{sample}/umi/{sample}"+EXT
		log:OUT+"/{sample}/umi/{sample}.log"
		conda:
			"../envs/umi_and_trim.yaml"
		params:
			trimfile=OUT+"/{sample}/umi/{sample}_fastx"+EXT,
			trimnumber=config["umi_tools"]["trimm"]["trim_supp"]+1,
			extra= umi_extract_extra,
			bc_pattern=config["umi_tools"]["extract"]["bc-pattern"],
			extract_method=config["umi_tools"]["extract"]["extract-method"]
		shell:
			"umi_tools extract --bc-pattern={params.bc_pattern} --extract-method={params.extract_method} --log={log} --stdout={params.trimfile} --stdin={input.fastq} {params.extra} && "
			"fastx_trimmer -f {params.trimnumber} -i {params.trimfile} -o {output.fastq} && "
			"rm {params.trimfile}"
	
else:
	rule umi_extract:
		input:
			fastq=IN+"/{sample}"+EXT
		output:
			fastq=OUT+"/{sample}/umi/{sample}"+EXT
		log:OUT+"/{sample}/umi/{sample}.log"
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
		bam=OUT+"/{prefix}/{genome}/bam/sorted/{prefix}.sorted.bam",
		bai=OUT+"/{prefix}/{genome}/bam/sorted/{prefix}.sorted.bai"
	output:
		bam=OUT+"/{prefix}/{genome}/bam/umidedup/{prefix}.umidedup.bam"
	conda:
		"../envs/umi_tools.yaml"
	params:
		stats=OUT+"/{prefix}/{genome}/bam/umidedup/{prefix}"
	shell:
		"umi_tools dedup -I {input.bam} --output-stats={params.stats} -S {output.bam}"
