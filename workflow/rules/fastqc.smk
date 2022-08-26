rule fastqc_raw:
	input:
		fastq=OUT+"/{prefix}/demultiplex/{prefix}"+config["fastq"]
	output:
		OUT+"/{prefix}/QC/RAW/{prefix}_raw_fastqc.html"
	params:
		dir=OUT+"/{prefix}/QC/RAW/",
		before=OUT+"/{prefix}/QC/RAW/{prefix}_fastqc.zip",
		after=OUT+"/{prefix}/QC/RAW/{prefix}_raw_fastqc.zip"
	threads: config["threads"]
	benchmark :
		OUT+"/{prefix}/benchmarks/fastqc_raw/{prefix}.txt"
	priority: 100
	conda: "../envs/fastqc.yaml"
	message : "##RUNNING : fastqc for {input.fastq}"
	shell: 
		"fastqc -q -t {threads} --outdir {params.dir} {input.fastq} && mv {before} {after}"



rule fastqc_bam:
	input:
		bam=OUT+"/{prefix}/mapping/bam/{samtype}/{prefix}.{samtype}.bam"
	output:
		OUT+"/{prefix}/QC/bam/{prefix}.{samtype}_fastqc.html"
	params:
		dir=OUT+"/{prefix}/QC/bam/"
	threads: config["threads"]
	benchmark :
		OUT+"/{prefix}/benchmarks/fastqc_bam/{prefix}.{samtype}.txt"
	priority: 10
	message : "##RUNNING : fastqc for {input.bam}"
	conda: "../envs/fastqc.yaml"
	shell: "fastqc -q -t {threads} --outdir {params.dir} {input.bam}"


rule multi_qc:
	input:
		expand(OUT+"/{sample}/QC/RAW/{sample}_raw_fastqc.zip",sample=config["samples"]),
		expand(OUT+"/{sample}/QC/bam/{sample}.{samtype}_fastqc.html",samtype=samtype,sample=config["samples"]),
		expand(OUT+"/{sample}/QC/STATS/{type}/{sample}.{samtype}.{type}",samtype=samtype,type = ["stats","idxstats","flagstat"],sample=config["samples"])
	output : 
		OUT+"/{sample}/QC/MULTIQC/"+OUT+"_multiqc_report.html"
	params:
		title = OUT,
		conf = config["fastqc"]["multi_qc_path"],
		output = OUT+"/{sample}/QC/MULTIQC",
		filename = OUT+"_multiqc_report.html"
	priority: 50
	message : "##RUNNING : MultiQC"
	conda: "../envs/fastqc.yaml"
	shell:
		"export LC_ALL=C.UTF-8 && "
		"export LANG=C.UTF-8 && "
		"multiqc {params.title} --config {params.conf} --title {params.title} -o {params.output} --filename {params.filename}"