suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(tidyverse)))
dataset_raw <- snakemake@input[[1]]
message(dataset_raw)

myBam <- scanBam(dataset_raw)[[1]]

myBam_seq_name <- tibble(names = myBam[["qname"]],seq = as.character(myBam[["seq"]]))

myBam_tibble <- GAlignments(seqnames =myBam[["rname"]],pos = myBam[["pos"]],cigar = myBam[["cigar"]],
                        strand = myBam[["strand"]],names = myBam[["qname"]])%>% as.data.frame() %>% tibble::rownames_to_column("names")  %>% as_tibble()


myBam_tibble <- myBam_tibble %>% left_join(myBam_seq_name)
myBam_tibble <- myBam_tibble %>% mutate(UMI = str_extract(names,"[ACGT]+$"))

stats <- myBam_tibble %>% dplyr::count(start,end,strand,UMI,seq) 
res <- myBam_tibble %>% dplyr::group_by(start,end,strand,UMI,seq) %>% slice(1) %>% pull(names)

res %>% enframe() %>% dplyr::select(value) %>% write_tsv(snakemake@output[[1]],col_names = F)
stats %>% write_tsv(snakemake@output[[2]])

