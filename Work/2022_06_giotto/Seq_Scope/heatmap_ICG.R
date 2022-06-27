heatmap_ICG <- function(df_ICG){
  
  uniq_ct <- union(df_ICG$receiver_cell_type, df_ICG$sender_cell_type)
  ct_num <- length(uniq_ct)
  
  # Up regulation
  df_ICG_up <- df_ICG[df_ICG$log2fc >= 0, ]
  
  icg_mat_up <- matrix(0, nrow = ct_num, ncol = ct_num)
  rownames(icg_mat_up) <- uniq_ct
  colnames(icg_mat_up) <- uniq_ct
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
      
      r_flag <- df_ICG_up$receiver_cell_type %in% uniq_ct[r_ct]
      s_flag <- df_ICG_up$sender_cell_type %in% uniq_ct[s_ct]
      
      t_flag <- r_flag & s_flag

      df_ICG_up_t <- df_ICG_up[t_flag, ]
      
      icg_mat_up[r_ct, s_ct] <- nrow(df_ICG_up_t)

    }
  }
  
  write.csv(icg_mat_up, file = "icg_mat_up.csv")
  
  pdf("ICGs_up.pdf")
  NMF::aheatmap(icg_mat_up, Rowv = NA, Colv = NA, color = "Reds", txt = icg_mat_up,
                main = "Giotto up-regulation ICGs")
  dev.off()
  
  # Down regulation
  df_ICG_down <- df_ICG[df_ICG$log2fc < 0, ]
  
  icg_mat_down <- matrix(0, nrow = ct_num, ncol = ct_num)
  rownames(icg_mat_down) <- uniq_ct
  colnames(icg_mat_down) <- uniq_ct
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
      
      r_flag <- df_ICG_down$receiver_cell_type %in% uniq_ct[r_ct]
      s_flag <- df_ICG_down$sender_cell_type %in% uniq_ct[s_ct]
      
      t_flag <- r_flag & s_flag
      
      df_ICG_down_t <- df_ICG_down[t_flag, ]
      
      icg_mat_down[r_ct, s_ct] <- nrow(df_ICG_down_t)
      
    }
  }
 
  write.csv(icg_mat_down, file = "icg_mat_down.csv")
  
  pdf("ICGs_down.pdf")
  NMF::aheatmap(icg_mat_down, Rowv = NA, Colv = NA, color = "Blues", txt = icg_mat_down,
                main = "Giotto down-regulation ICGs")
  dev.off()
  
}