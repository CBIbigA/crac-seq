suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(tidyverse)))
dataset_raw <- snakemake@input[[1]]
message(dataset_raw)

BF <- BamFile(dataset_raw) ;
param <- ScanBamParam(what="seq")


myBam_tibble <- readGAlignments(BF,use.names = T, param=param) %>% as.data.frame() %>% tibble::rownames_to_column("names")  %>% as_tibble()

myBam_tibble <- myBam_tibble %>% mutate(UMI = str_extract(names,"[ACGT]+$"))

stats <- myBam_tibble %>% dplyr::count(start,end,strand,UMI,seq) 
res <- myBam_tibble %>% dplyr::group_by(start,end,strand,UMI,seq) %>% slice(1) %>% pull(names)

res %>% enframe() %>% dplyr::select(value) %>% write_tsv(snakemake@output[[1]],col_names = F)
stats %>% write_tsv(snakemake@output[[2]])

