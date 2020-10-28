# Get alignment statistics for Hi-C data

## load packages
library("readr")
library("magrittr")
library("Rsamtools")
library("furrr")
library("R.utils")

## get metadata
metadata <- read_tsv(file.path("metadata.txt"), col_names = TRUE) %>%
  tidyr::gather(key = "file_type", value = "file_name", "run1_read1":"run4_read2") %>%
  tidyr::extract(file_type, into = c("Technical replicate", "Read"), regex = "run([[:alnum:]]+)_read([[:alnum:]]+)") %>%
  dplyr::filter(!is.na(file_name))

metadata <- metadata %>%
  dplyr::rename(sample_name = sample)

metadata <- metadata %>%
  dplyr::mutate(fastq_location = file.path("data", "fastq", file_name),
                bam_location = file.path("data", "hic", sample_name, "sam", 
                                         gsub(".fastq.gz|.fq.gz", "_sort.bam", file_name)))

## set up parameters for reading in bam file data
scanbam_param1 <- Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
                                          mapqFilter = 3, tag="XA")

countBam2 <- function(file, param = ScanBamParam()){
  bam <- scanBam(file, param = param)
  return(sum(is.na(bam[[1]]$tag$XA)))
}

countLines2 <- function(file){
  countLines(file)[1]
}

## calculate statistics
plan(multiprocess)
sample_stats <- metadata %>%
  dplyr::mutate(fastq_lines = furrr::future_map_int(fastq_location, countLines2),
                fastq_records = fastq_lines / 4,
                bam_records = furrr::future_map_int(bam_location, countBam2, param = scanbam_param1))

## write stats to file
write_tsv(sample_stats, path = file.path("data", "sample_alignment_stats_bwa.txt"))

sessionInfo()
