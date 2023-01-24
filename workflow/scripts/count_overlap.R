suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(plyranges)))


monbam <- snakemake@input[["bam"]]
monbai <- snakemake@input[["bai"]]

BF <- BamFile(monbam,index=monbai) ;
## Read file
GA <- readGAlignments(BF) ;
GR <- GRanges(GA)


monbed <- read_bed(snakemake@input[["bed"]])

monbed$score <- plyranges::count_overlaps(monbed,GR)

export.bedGraph(monbed,snakemake@output[[1]])