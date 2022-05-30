rule novoindex:
	input:
		genome="{0}/{1}{2}".format(config["genome"]["dir"],config["genome"]["name"],config["genome"]["ext"])
	output:
		index="{0}/{1}.index".format(config["genome"]["dir"],config["genome"]["name"])
	threads:config["threads"]
	conda:
		"../envs/novoalign.yaml"
	shell:
		"novoindex -t {threads} {output.index} {input.genome}"


rule novoalign:
	input:
		index="{0}/{1}.index".format(config["genome"]["dir"],config["genome"]["name"]),
		fastq=OUT+"/demultiplex/{prefix}"+in_ext
	threads:1
	output:
		OUT+"/mapping/sam/{prefix}.sam"
	conda:
		"../envs/novoalign.yaml"
	shell:
		"novoalign -c {threads} -f {input.fastq} -d {input.index} -r Random -o SAM > {output}"