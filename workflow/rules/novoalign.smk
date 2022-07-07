

rule novoindex:
	input:
		genome=getGenome
	output:
		index="OUT/{genome}/{genome}.index"
	threads:config["threads"]
	conda:
		"../envs/novoalign.yaml"
	shell:
		"novoindex -t {threads} {output.index} {input.genome}"



if not config["cutadapt"]["adapters"]:
	in_ext = config["fastq"]
else:
	in_ext = ".trimmed"+config["fastq"]


def conditionalInputNovoAlign(wildcards):
	if config["demultiplexing"]["demultiplex"]:
		return(OUT+"/demultiplex/"+wildcards.sample+in_ext)
	else:
		return(OUT+"/umi/"+wildcards.sample+in_ext)

rule novoalign:
	input:
		genome=getGenome,
		index="OUT/{genome}/{genome}.index",
		fastq=conditionalInputNovoAlign
	threads:1
	output:
		OUT+"/{genome}/mapping/sam/{sample}.sam"
	conda:
		"../envs/novoalign.yaml"
	shell:
		"novoalign -c {threads} -f {input.fastq} -d {input.index} -r Random -o SAM > {output}"