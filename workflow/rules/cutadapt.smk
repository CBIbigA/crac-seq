def getOutFormat(wildcards):
  return(OUT+"/demultiplex/"+wildcards.sample+"_{name}"+config["fastq"])



rule cutadapt_demultiplex:
	input:
		barcode = OUT+"/barcode.fa",
		fastq = OUT+"/umi/{sample}"+config["fastq"]
	output:
		expand(OUT+"/demultiplex/{{sample}}_{barcode}"+config["fastq"],barcode=config["demultiplexing"]["barcode"])
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
		fastq=OUT+"/{type}/{prefix}"+config["fastq"]
	output:
		fastq=OUT+"/{type}/{prefix}.trimmed"+config["fastq"],
		qc=OUT+"/{type}/{prefix}.trimmed.qc"
	threads:config["threads"]
	params:
		adapters = config["cutadapt"]["adapters"],
		extra = config["cutadapt"]["extra"]
	conda:
		"../envs/cutadapt.yaml"
	shell:
		"cutadapt {params.adapters} {params.extra} -j {threads} -o {output.fastq} {input} > {output.qc}"
