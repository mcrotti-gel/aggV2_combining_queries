#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input_table <- args[1]

par_var_tbl <- read_tsv(input_table, col_names = FALSE)

par_tbl <- par_var_tbl %>%
 group_by(X1) %>%
 count() %>%
 arrange(desc(n)) %>%
 rename(platekey = X1, n_variants = n)

var_tbl <- par_var_tbl %>%
 mutate(variant = paste(X2,X3,X4,X5, sep = "_")) %>%
 group_by(variant) %>%
 count() %>%
 arrange(desc(n)) %>%
 rename(n_platekey = n)

out_file_name_platekey <- paste0(sub("_results.tsv","",basename(input_table)),"_platekey_summary.tsv")
out_file_name_variant <- paste0(sub("_results.tsv","",basename(input_table)),"_variant_summary.tsv")

write.table(par_tbl, file = out_file_name_platekey, quote = F, sep = "\t", row.names = F)
write.table(var_tbl, file = out_file_name_variant, quote = F, sep = "\t", row.names = F)