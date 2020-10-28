library("readr")
library("dplyr")

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
threshold <- args[2]

message("Working on ", input_file)
message("Using boundary score threshold: ", threshold)

boundaries_df <- read_tsv(input_file, col_names = c("chr", "start", "end", "name", "score", "strand")) %>% 
  unique() %>%
  filter(!is.na(score))

total <- nrow(boundaries_df)
pass_threshold <- sum(boundaries_df$score > threshold)
percent_pass <- signif(pass_threshold/total*100, 3)

message("There are ", total, " unique boundaries")
message(pass_threshold, " boundaries pass the threshold (", percent_pass, "%)")

boundaries_filtered <- boundaries_df %>%
  unique() %>%
  filter(!is.na(score)) %>%
  filter(score > threshold)

output_file <- gsub(".bed", paste0("_filtered", threshold, ".bed"),
                    input_file)
message("Writing filtered boundaries to output file ", output_file)

boundaries_filtered %>%
  select(chr, start, end, name, score) %>% 
  write_tsv(path = output_file, col_names = FALSE)

sessionInfo()


