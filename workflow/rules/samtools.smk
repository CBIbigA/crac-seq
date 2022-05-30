rule samtools_faidx:
	input:
		genome="{0}/{1}{2}".format(config["genome"]["dir"],config["genome"]["name"],config["genome"]["ext"])
	output:
		"{0}/{1}{2}.fai".format(config["genome"]["dir"],config["genome"]["name"],config["genome"]["ext"])
	conda: "../envs/samtools.yaml"
	shell:
		"samtools faidx -o {output} {input}"

rule sam_to_bam:
	input:
		sam=OUT+"/mapping/sam/{prefix}.sam",
		genome="{0}/{1}{2}.fai".format(config["genome"]["dir"],config["genome"]["name"],config["genome"]["ext"])
	output:
		OUT+"/mapping/bam/raw/{prefix}.bam"
	params:
		quality=config["samtools"]["quality"],
		custom=config["samtools"]["custom"]
	benchmark :
		OUT+"/benchmarks/sam_to_bam/{prefix}.txt"
	priority: 50
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	message: "##RUNNING : samtools view for {input.sam}"
	shell:
		"samtools view "
		"{params.custom} -@ {threads} "
		"-b -S "
		"-q {params.quality} "
		"-t {input.genome} "
		"-o {output} "
		"{input.sam}"



rule samtools_sortn:
	input:
		bam=OUT+"/mapping/bam/raw/{prefix}.bam"
	output:
		temp(OUT+"/mapping/bam/sorted/{prefix}.nsorted.bam")
	benchmark :
		OUT+"/benchmarks/samtools_sortn/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools sort -n {input.bam}"
	conda: "../envs/samtools.yaml"
	threads: config["threads"]
	shell:
		"samtools sort -@ {threads} -n -o {output} "
		"{input.bam} "

rule samtools_fixmate:
	input:
		bam=OUT+"/mapping/bam/sorted/{prefix}.nsorted.bam"
	output:
		temp(OUT+"/mapping/bam/sorted/{prefix}.fixmate.bam")
	benchmark :
		OUT+"/benchmarks/samtools_fixmate/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools fixmate -m {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools fixmate -m -@ {threads} "
		"{input.bam} "
		"{output}"

rule samtools_sort:
	input:
		bam=OUT+"/mapping/bam/sorted/{prefix}.fixmate.bam"
	output:
		OUT+"/mapping/bam/sorted/{prefix}.sorted.bam"
	benchmark :
		OUT+"/benchmarks/samtools_sort/{prefix}.txt"
	priority: 50
	threads: config["threads"]
	message: "##RUNNING : samtools sort {input.bam}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools sort -@ {threads} "
		"-o {output} "
		"{input.bam}"


rule samtools_markdups:
	input:
		bam=OUT+"/mapping/bam/sorted/{prefix}.sorted.bam"
	output:
		OUT+"/mapping/bam/markdup/{prefix}.markdup.bam"
	benchmark :
		OUT+"/benchmarks/samtools_markdups/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools markdup {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools markdup -@ {threads} "
		"{input.bam} "
		"{output}"

rule samtools_rmdups:
	input:
		bam=OUT+"/mapping/bam/markdup/{prefix}.markdup.bam",
		bai=OUT+"/mapping/bam/markdup/{prefix}.markdup.bai"
	output:
		OUT+"/mapping/bam/rmdups/{prefix}.rmdups.bam"
	benchmark :
		OUT+"/benchmarks/samtools_rmdups/{prefix}.txt"
	priority: 50
	message: "##RUNNING : samtools view -F 1024 {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools view -@ {threads} -F 1024 -b "
		"{input.bam} "
		"> {output}"


rule samtools_index:
	input:
		bam=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam"
	output:
		OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	benchmark :
		OUT+"/benchmarks/samtools_index/{prefix}.{samtype}.txt"
	priority: 50
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	message: "##RUNNING : samtools index {input}"
	shell:
		"samtools index -@ {threads} {input} {output}"


rule samtools_stats:
	input:
		bam=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/QC/STATS/stats/{prefix}.{samtype}.stats"
	benchmark :
		OUT+"/benchmarks/samtools_stats/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools stats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools stats {input.bam} > {output}"

rule samtools_idxstats:
	input:
		bam=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/QC/STATS/idxstats/{prefix}.{samtype}.idxstats"
	benchmark :
		OUT+"/benchmarks/samtools_idxstats/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools idxstats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools idxstats {input.bam} > {output}"

rule samtools_flagstat:
	input:
		bam=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/QC/STATS/flagstat/{prefix}.{samtype}.flagstat"
	benchmark :
		OUT+"/benchmarks/samtools_flagstat/{prefix}.{samtype}.txt"
	priority: 50
	message: "##RUNNING : samtools flagstat {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools flagstat {input.bam} > {output}"
