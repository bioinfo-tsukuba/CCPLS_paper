collectICG <- function(tsv_prefix, findICG_result){
  
  sig_gene_df <- findICG_result$CPGscores
  
  receiver_cell_type <- sig_gene_df$cell_type
  sender_cell_type <- sig_gene_df$int_cell_type
  sig_gene_symbol <- sig_gene_df$genes
  log2fc <- sig_gene_df$log2fc
  
  df_save_1 <- as_tibble(sig_gene_symbol)
  names(df_save_1) <- "sig_gene_symbol"
  
  df_save_1 %>% dplyr::mutate(receiver_cell_type = receiver_cell_type) %>%
                      dplyr::mutate(sender_cell_type = sender_cell_type) %>%
                      dplyr::mutate(log2fc = log2fc) -> df_save_2
  
  write.table(df_save_2, paste0(tsv_prefix, "_ICG.tsv"), sep = "\t", col.names = FALSE)

  return(df_save_2)
  
}