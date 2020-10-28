gene_ids_key <- rtracklayer::import.gff("../../external_data/flybase/dmel-all-r6.30.gtf.gz") %>% 
  as.data.frame() %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_symbol)

txdb <- makeTxDbFromGFF("../../external_data/flybase/dmel-all-r6.30.gtf.gz")
transcripts_gr <- transcripts(txdb, columns = c("TXNAME", "GENEID"))
mcols(transcripts_gr)$GENEID <- unlist(mcols(transcripts_gr)$GENEID)

mcols(transcripts_gr) <- left_join(as.data.frame(mcols(transcripts_gr)), 
                                   gene_ids_key, by = c("GENEID" = "gene_id"))

transcripts_gr <- transcripts_gr[seqnames(transcripts_gr) %in% c("2L", "2R", "3L", "3R", "4", "X", "Y")]

# potential transcripts to exclude
# sapply(transcripts_gr$gene_symbol[grepl(":", transcripts_gr$gene_symbol)], function(str){
#   strsplit(str, ":", fixed = TRUE)[[1]][1]
# }) %>%  table()

tx_types_to_exclude <- c("tRNA:", "rRNA:", "rRNA-Psi:", "asRNA:", "snoRNA:", "scaRNA:", "snRNA:")

transcripts_gr <- transcripts_gr[!grepl(pattern = paste(tx_types_to_exclude, collapse = "|"), 
                                        transcripts_gr$gene_symbol)]

promoters_gr <- promoters(transcripts_gr, upstream = 1, downstream = 0)


promoters_gr[promoters_gr$gene_symbol %in% c("twi", "sna", "if", "NetA", "sog", "Doc1", "pnr", "C15")]

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "twi"][3], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "sna"], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "if"][3], width = 2000, fix = "center"))

as.character(resize(range(promoters_gr[promoters_gr$gene_symbol == "NetA"]), width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "sog"][1], width = 2000, fix = "center"))

as.character(resize(range(promoters_gr[promoters_gr$gene_symbol == "Doc1"]), width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "pnr"][3], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "C15"], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "brk"], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "sog"][3], width = 2000, fix = "center"))


as.character(resize(promoters_gr[promoters_gr$gene_symbol == "ap"][2], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "pdm2"][4], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "Con"][2], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "eya"][3], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "stumps"][2], width = 2000, fix = "center"))

as.character(resize(sort(promoters_gr[promoters_gr$gene_symbol == "Mef2"])[6], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "sli"][1], width = 2000, fix = "center"))

as.character(resize(promoters_gr[promoters_gr$gene_symbol == "slp1"], width = 2000, fix = "center"))

as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "Abd-B"]))[5], width = 2000, fix = "center"))

as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "E2f1"]))[5], width = 2000, fix = "center"))

## housekeeping genes!

as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "RpS12"])), width = 2000, fix = "center"))
as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "eEF1delta"])), width = 2000, fix = "center"))
as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "x16"])), width = 2000, fix = "center"))
as.character(resize(sort(unique(promoters_gr[promoters_gr$gene_symbol == "Nipped-B"])), width = 2000, fix = "center"))
