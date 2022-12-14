get_coverage <- function(my.GR){
  CovBam <- coverage(my.GR) ;
  # nCovBam <- CovBam * nf ;
  # seqlengths( nCovBam ) <- seqlens[names(seqlengths( nCovBam ) )] ;
  res_coverage <- split(mes_pos_genome,seqnames(mes_pos_genome)) %>% lapply(function(x){
    x$score <- Views(CovBam[[seqnames(x[1])@values]],start=start(x),end=end(x))  %>% subject() %>% as.numeric() %>% as.numeric()
    return(x)
  }) %>% GRangesList() %>% unlist()
  return(res_coverage)
}

#.libPaths(c("/usr/local/lib/R/site-library",.libPaths()))
library(Rsamtools) ; # ScanBam
library(GenomicAlignments) ; # GAlignement
library( plyranges ) ; #Export bw
library(rtracklayer)
library( Biostrings ) ; #Export bw
library(stringr)
library(dplyr)
seqlens <- snakemake@input[["genome"]]


filename <- snakemake@input[["bam"]];

##Normalisation 
BF <- BamFile(filename) ;
nbReads <- countBam(BF)$records ;
nf <- as.numeric(1000000/nbReads) ;

##Seqlens
seqlens <- read.table(seqlens,sep="\t",header = F)
seqlens <- setNames(seqlens[,2], seqlens[,1])
mes_pos_genome <- tileGenome(seqlens,tilewidth = 1,cut.last.tile.in.chrom = T)


## Read file
GA <- readGAlignments(BF) ;
GR <- GRanges(GA)
seqlengths( GR ) <- seqlens[names(seqlengths( GR ) )] ;
##Normalized bedGraph
# CovBam <- coverage(GA) ;
# nCovBam <- CovBam * nf ;
# seqlengths( nCovBam ) <- seqlens[names(seqlengths( nCovBam ) )] ;
# export.bedGraph(nCovBam,snakemake@output["normalized"]) ;
# export.bedGraph(CovBam,snakemake@output["unnormalized"]) ;




message("Export bedGRaph")
GR %>% get_coverage() %>% export.bedGraph(snakemake@output[[1]])
##3 prime position
message("Export bedGRaph on 3prime position only")
GR %>% anchor_3p() %>% mutate(width=1) %>% get_coverage() %>%  export.bedGraph(snakemake@output[[2]])


##Take only substitution reads

message("Export bedGRaph on 3prime position only for read with substitution")
GR.sub <- GA %>% as.data.frame() %>% 
  mutate(Substitution_prime3 = case_when(
    strand =="-" & str_detect(cigar,"^[0-9]+S.+") ~ TRUE,
    strand =="+" & str_detect(cigar,".+[0-9]+S$") ~ TRUE,
    TRUE ~ FALSE
  )) %>%  filter(Substitution_prime3) %>% as_granges()
seqlengths( GR.sub ) <- seqlens[names(seqlengths( GR.sub ) )] ;
GR.sub %>% anchor_3p() %>% mutate(width=1) %>% get_coverage() %>% export.bedGraph(snakemake@output[[3]])

message("Export bedGRaph on 3prime position only for read with substitution")
GR.sub <- GA %>% as.data.frame() %>% 
  mutate(Substitution_prime3 = case_when(
    strand =="+" & str_detect(cigar,"^[0-9]+S.+") ~ TRUE,
    strand =="-" & str_detect(cigar,".+[0-9]+S$") ~ TRUE,
    TRUE ~ FALSE
  )) %>%  filter(Substitution_prime3) %>% as_granges()
seqlengths( GR.sub ) <- seqlens[names(seqlengths( GR.sub ) )] ;
GR.sub %>% anchor_3p() %>% mutate(width=1) %>% get_coverage() %>% export.bedGraph(snakemake@output[[4]])


