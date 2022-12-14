rule R_umi_dedup:
  input:
    bam=OUT+"/{prefix}/{genome}/bam/sorted/{prefix}.sorted.bam",
    bai=OUT+"/{prefix}/{genome}/bam/sorted/{prefix}.sorted.bai"
  output:
    OUT+"/{prefix}/{genome}/bam/Rumidedup/{prefix}_list.txt",
    OUT+"/{prefix}/{genome}/bam/Rumidedup/{prefix}_stats.tsv"
  threads:4
  conda:
    "../envs/Rumidedup.yaml"
  script:
    "../scripts/RemoveDuplicates_UMI.R"


rule picard_filter_reads:
  input:
    bam=OUT+"/{prefix}/{genome}/bam/sorted/{prefix}.sorted.bam",
    listreads=OUT+"/{prefix}/{genome}/bam/Rumidedup/{prefix}_list.txt"
  output:
    bam=OUT+"/{prefix}/{genome}/bam/Rumidedup/{prefix}.Rumidedup.bam"
  conda:
    "../envs/picard.yaml"
  shell:
    "picard FilterSamReads I={input.bam} O={output.bam} READ_LIST_FILE={input.listreads} FILTER=includeReadList"


