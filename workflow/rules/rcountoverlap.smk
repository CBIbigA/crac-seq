rule generate_coverage:
  input:
    bam=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bam",
    bai=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bai",
    genome=OUT+"/{genome}_chrom.size"
  output:
    OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.bedGraph",
    OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.3primebedGraph",
    OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.3prime_3prime_substitutionbedGraph",
    OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.3prime_5prime_substitutionbedGraph"
  conda: "../envs/coverage.yaml"
  script:
    "../scripts/generate_coverage_perBase.R"

def getbed(wildcards):
  return(config["genome"][wildcards.genome]["bed"])


# rule generate_coverage_for_rDNA:
#   input:
#     bam=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bam",
#     bai=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bai",
#     bed=getbed
#   output:
#     OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.bedReporting.tsv"
#   conda: "../envs/bedtools.yaml"
#   #threads: 4
#   shell:
#     "bedtools coverage -nonamecheck -a {input.bed} -b {input.bam} > {output}"

rule R_count_overlap:
  input:
    bam=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bam",
    bai=OUT+"/{prefix}/{genome}/bam/{samtype}/{prefix}.{samtype}.bai",
    bed=getbed
  output:
    OUT+"/{prefix}/{genome}/coverage/{samtype}/{prefix}.{samtype}_{genome}.bedReporting.tsv"
  threads:4
  conda:
    "../envs/Rumidedup.yaml"
  script:
    "../scripts/count_overlap.R"
