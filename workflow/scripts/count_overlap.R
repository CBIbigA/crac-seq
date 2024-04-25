suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(plyranges)))
suppressWarnings(suppressMessages(library(rtracklayer)))


monbam <- snakemake@input[["bam"]]
monbai <- snakemake@input[["bai"]]

BF <- BamFile(monbam,index=monbai) ;
## Read file
GA <- readGAlignments(BF) ;
GR <- GRanges(GA)


monbed <- read_bed(snakemake@input[["bed"]])

monbed$score <- plyranges::count_overlaps(monbed,GR)

#export.bedGraph(monbed,snakemake@output[[1]])
write.table( as.data.frame( monbed )[, c(1, 2, 3, 7, 6 )], file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
