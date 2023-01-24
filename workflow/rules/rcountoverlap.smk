def getbed(wildcards):
  return(config["genome"][wildcards.genome]["bed"])

rule R_count_overlap:
  input:
    bam=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bam",
    bai=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bai",
    bed=getbed
  output:
    OUT+"/{prefix}/{genome}/Rcoverage/{samtype}/{prefix}.{samtype}_{genome}.RbedReporting.tsv"
  threads:4
  conda:
    "../envs/Rumidedup.yaml"
  script:
    "../scripts/count_overlap.R"