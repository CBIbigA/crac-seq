rule samtools_faidx:
	input:
		genome=getGenome
	output:
		"OUT/{genome}/{genome}.fai"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools faidx -o {output} {input}"


rule makeChromSize:
	input:
		genome=getGenome
	output:
		chromsize=OUT+"/{genome}_chrom.size"
	conda: "../envs/pipeline.yaml"
	shell:
		"faidx {input} -i chromsizes > {output}"

rule sam_to_bam:
	input:
		sam=OUT+"/{genome}/mapping/sam/{prefix}.sam",
		genome="OUT/{genome}/{genome}.fai"
	output:
		OUT+"/{genome}/mapping/bam/raw/{prefix}.bam"
	params:
		quality=config["samtools"]["quality"],
		custom=config["samtools"]["custom"]
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
		bam=OUT+"/{genome}/mapping/bam/raw/{prefix}.bam"
	output:
		OUT+"/{genome}/mapping/bam/sorted/{prefix}.nsorted.bam"
	message: "##RUNNING : samtools sort -n {input.bam}"
	conda: "../envs/samtools.yaml"
	threads: config["threads"]
	shell:
		"samtools sort -@ {threads} -n -o {output} "
		"{input.bam} "

rule samtools_fixmate:
	input:
		bam=OUT+"/{genome}/mapping/bam/sorted/{prefix}.nsorted.bam"
	output:
		OUT+"/{genome}/mapping/bam/sorted/{prefix}.fixmate.bam"
	message: "##RUNNING : samtools fixmate -m {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools fixmate -m -@ {threads} "
		"{input.bam} "
		"{output}"

rule samtools_sort:
	input:
		bam=OUT+"/{genome}/mapping/bam/sorted/{prefix}.fixmate.bam"
	output:
		OUT+"/{genome}/mapping/bam/sorted/{prefix}.sorted.bam"
	threads: config["threads"]
	message: "##RUNNING : samtools sort {input.bam}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools sort -@ {threads} "
		"-o {output} "
		"{input.bam}"


rule samtools_markdups:
	input:
		bam=OUT+"/{genome}/mapping/bam/sorted/{prefix}.sorted.bam"
	output:
		OUT+"/{genome}/mapping/bam/markdup/{prefix}.markdup.bam"
	message: "##RUNNING : samtools markdup {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools markdup -@ {threads} "
		"{input.bam} "
		"{output}"

rule samtools_rmdups:
	input:
		bam=OUT+"/{genome}/mapping/bam/markdup/{prefix}.markdup.bam",
		bai=OUT+"/{genome}/mapping/bam/markdup/{prefix}.markdup.bai"
	output:
		OUT+"/{genome}/mapping/bam/rmdups/{prefix}.rmdups.bam"
	message: "##RUNNING : samtools view -F 1024 {input.bam}"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	shell:
		"samtools view -@ {threads} -F 1024 -b "
		"{input.bam} "
		"> {output}"


rule samtools_index:
	input:
		bam=OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bam"
	output:
		OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	threads: config["threads"]
	conda: "../envs/samtools.yaml"
	message: "##RUNNING : samtools index {input}"
	shell:
		"samtools index -@ {threads} {input} {output}"


rule samtools_stats:
	input:
		bam=OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/{genome}/QC/STATS/stats/{prefix}.{samtype}.stats"
	message: "##RUNNING : samtools stats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools stats {input.bam} > {output}"

rule samtools_idxstats:
	input:
		bam=OUT+"//{genome}mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/{genome}/QC/STATS/idxstats/{prefix}.{samtype}.idxstats"
	message: "##RUNNING : samtools idxstats {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools idxstats {input.bam} > {output}"

rule samtools_flagstat:
	input:
		bam=OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bam",
		bai=OUT+"/{genome}/mapping/bam/{samtype}/{prefix}.{samtype}.bai"
	output:
		OUT+"/{genome}/QC/STATS/flagstat/{prefix}.{samtype}.flagstat"
	message: "##RUNNING : samtools flagstat {input}"
	conda: "../envs/samtools.yaml"
	shell:
		"samtools flagstat {input.bam} > {output}"
