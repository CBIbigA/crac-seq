def getOutFormat(wildcards):
  return(OUT+"/"+wildcards.sample+"_{name}"+"/demultiplex/"+wildcards.sample+"_{name}"+config["fastq"])



rule cutadapt_demultiplex:
	input:
		barcode = OUT+"/barcode.fa",
		fastq = OUT+"/{sample}/umi/{sample}"+config["fastq"]
	output:
		expand(OUT+"/{{sample}}_{barcode}/demultiplex/{{sample}}_{barcode}"+config["fastq"],barcode=config["demultiplexing"]["barcode"]),
		OUT+"/{sample}_unknown/demultiplex/{sample}_unknown"+config["fastq"]
	log:OUT+"/demultiplex/{sample}.demultiplex.qc"
	threads:config["threads"]
	params:
		extra = config["demultiplexing"]["extra"],
		outformat= getOutFormat,
		errorate = config["demultiplexing"]["errorate"]
	conda:
		"../envs/cutadapt.yaml"
	shell:
		"cutadapt -e {params.errorate}  -g file:{input.barcode} -o '{params.outformat}' {input.fastq} > {log}"


rule cutadapt:
	input:
		fastq=OUT+"/{prefix}/{type}/{prefix}"+config["fastq"]
	output:
		fastq=OUT+"/{prefix}/{type}/{prefix}.trimmed"+config["fastq"],
		qc=OUT+"/{prefix}/{type}/{prefix}.trimmed.qc"
	threads:config["threads"]
	params:
		adapters = config["cutadapt"]["adapters"],
		extra = config["cutadapt"]["extra"]
	conda:
		"../envs/cutadapt.yaml"
	shell:
		"cutadapt {params.adapters} {params.extra} -j {threads} -o {output.fastq} {input} > {output.qc}"
