#.libPaths(c("/usr/local/lib/R/site-library",.libPaths()))
library(Rsamtools) ; # ScanBam
library(GenomicAlignments) ; # GAlignement
library( plyranges ) ; #Export bw
library(rtracklayer)
library( Biostrings ) ; #Export bw

seqlens <- snakemake@input[["genome"]]
filename <- snakemake@input[["bam"]];
##Normalisation 
BF <- BamFile(filename) ;
nbReads <- countBam(BF)$records ;
nf <- as.numeric(1000000/nbReads) ;

##Seqlens
seqlens <- read.table(seqlens,sep="\t",header = F)
seqlens <- setNames(seqlens[,2], seqlens[,1])

## Read file
GA <- readGAlignments(BF) ;

##Normalized bedGraph
# CovBam <- coverage(GA) ;
# nCovBam <- CovBam * nf ;
# seqlengths( nCovBam ) <- seqlens[names(seqlengths( nCovBam ) )] ;
# export.bedGraph(nCovBam,snakemake@output["normalized"]) ;
# export.bedGraph(CovBam,snakemake@output["unnormalized"]) ;

##3 prime position

prime_3 <- GRanges(GA)
prime_3 <- mutate(anchor_3p(GRanges(GA)),width=1)

CovBam <- coverage(prime_3) ;
nCovBam <- CovBam * nf ;
seqlengths( nCovBam ) <- seqlens[names(seqlengths( nCovBam ) )] ;

# export.bedGraph(nCovBam,snakemake@output["3prime_normalized"]) ;
export.bedGraph(CovBam,snakemake@output[[1]])
