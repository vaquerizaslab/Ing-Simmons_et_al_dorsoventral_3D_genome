np_files <- list.files("./chipseq_peaks", pattern = "*.narrowPeak", full.names = TRUE)

lapply(np_files, function(f){
  peaks <- read.table(f, col.names = c("chr", "start", "end", "name", "score", "strand",
                              "fc", "pval", "qval", "summit"))
  peaks <- peaks[peaks$chr %in% c("2L", "2R", "3L", "3R", "4", "X", "Y"), ]
  outfile <- gsub(".narrowPeak", "_filtered.narrowPeak", fixed = TRUE, f)
  message(outfile)
  write.table(peaks, file = outfile, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t")
})


bp_files <- list.files("./chipseq_peaks", pattern = "*.broadPeak", full.names = TRUE)

lapply(bp_files, function(f){
  peaks <- read.table(f, col.names = c("chr", "start", "end", "name", "score", "strand",
                                       "fc", "pval", "qval"))
  peaks <- peaks[peaks$chr %in% c("2L", "2R", "3L", "3R", "4", "X", "Y"), ]
  outfile <- gsub(".broadPeak", "_filtered.broadPeak", fixed = TRUE, f)
  message(outfile)
  write.table(peaks, file = outfile, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t")
})
