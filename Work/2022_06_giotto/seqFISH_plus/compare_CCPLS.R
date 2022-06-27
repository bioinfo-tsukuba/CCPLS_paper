compare_CCPLS <- function(df_ICG, res.main){
  
  df_ICG$receiver_cell_type <- gsub(" ", ".", df_ICG$receiver_cell_type)
  df_ICG$sender_cell_type <- gsub(" ", ".", df_ICG$sender_cell_type)
  
  uniq_ct <- res.main$res.sep.mat$cell_type_list
  ct_num <- length(uniq_ct)
  
  ## Collect ICGs
  
  # Up
  ICG_symbol_list_up <- vector("list", ct_num * ct_num)
  df_ICG_up <- df_ICG[df_ICG$log2fc >= 0, ]  
  
  cnt <- 0
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
      
      cnt <- cnt + 1
      
      r_flag <- df_ICG_up$receiver_cell_type %in% uniq_ct[r_ct]
      s_flag <- df_ICG_up$sender_cell_type %in% uniq_ct[s_ct]
      t_flag <- r_flag & s_flag
      df_ICG_up_t <- df_ICG_up[t_flag, ]
      ICG_symbol_list_up[[cnt]] <- df_ICG_up_t$sig_gene_symbol

    }
  }
  
  # Down
  ICG_symbol_list_down <- vector("list", ct_num * ct_num)
  df_ICG_down <- df_ICG[df_ICG$log2fc < 0, ]  
  
  cnt <- 0
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
      
      cnt <- cnt + 1
      
      r_flag <- df_ICG_down$receiver_cell_type %in% uniq_ct[r_ct]
      s_flag <- df_ICG_down$sender_cell_type %in% uniq_ct[s_ct]
      t_flag <- r_flag & s_flag
      df_ICG_down_t <- df_ICG_down[t_flag, ]
      ICG_symbol_list_down[[cnt]] <- df_ICG_down_t$sig_gene_symbol
      
    }
  }
  
  ## Collect CCPLS
  
  # Up
  CCPLS_symbol_list_up <- vector("list", ct_num * ct_num)

  cnt <- 0
  
  for (r_ct in 1:ct_num){
    
    sig_mat <- res.main$res.sel.var$sig_coef_mat_list[[r_ct]]
    s_cand <- rownames(sig_mat)
    
    for (s_ct in 1:ct_num){
      
      cnt <- cnt + 1
      
      if (uniq_ct[s_ct] %in% s_cand){
        s_ct_row <- uniq_ct[s_ct]
        t_coef <- sig_mat[s_ct_row, ]
        t_symbol <- names(t_coef[t_coef > 0])
        CCPLS_symbol_list_up[[cnt]] <- t_symbol
      } else {
        CCPLS_symbol_list_up[[cnt]] <- c()
      }
      
    }
  }
  
  # Down
  CCPLS_symbol_list_down <- vector("list", ct_num * ct_num)
  
  cnt <- 0
  
  for (r_ct in 1:ct_num){
    
    sig_mat <- res.main$res.sel.var$sig_coef_mat_list[[r_ct]]
    s_cand <- rownames(sig_mat)
    
    for (s_ct in 1:ct_num){
      
      cnt <- cnt + 1
      
      if (uniq_ct[s_ct] %in% s_cand){
        s_ct_row <- uniq_ct[s_ct]
        t_coef <- sig_mat[s_ct_row, ]
        t_symbol <- names(t_coef[t_coef < 0])
        CCPLS_symbol_list_down[[cnt]] <- t_symbol
      } else {
        CCPLS_symbol_list_down[[cnt]] <- c()
      }
      
    }
  }
  
  ## Compare
  # Up
  ov_mat_up <- matrix(0, nrow = ct_num, ncol = ct_num)
  rownames(ov_mat_up) <- uniq_ct
  colnames(ov_mat_up) <- uniq_ct
  
  cnt <- 0
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
     
      cnt <- cnt + 1
      
      ICG <- ICG_symbol_list_up[[cnt]]
      CCPLS <- CCPLS_symbol_list_up[[cnt]]

      ov_mat_up[r_ct, s_ct] <- length(intersect(ICG, CCPLS))
       
    }
  }
  
  if (sum(ov_mat_up) > 0){
    pdf("ov_mat_up.pdf")
    NMF::aheatmap(ov_mat_up, Rowv = NA, Colv = NA, color = "Reds", txt = ov_mat_up,
                  main = "Overlap between Giotto findICG and CCPLS")
    dev.off()
  } else {
    print("No ovelapping in up-regulated genes.")
  }
  
  # Down
  ov_mat_down <- matrix(0, nrow = ct_num, ncol = ct_num)
  rownames(ov_mat_down) <- uniq_ct
  colnames(ov_mat_down) <- uniq_ct
  
  cnt <- 0
  
  for (r_ct in 1:ct_num){
    for (s_ct in 1:ct_num){
      
      cnt <- cnt + 1
      
      ICG <- ICG_symbol_list_down[[cnt]]
      CCPLS <- CCPLS_symbol_list_down[[cnt]]
      
      ov_mat_down[r_ct, s_ct] <- length(intersect(ICG, CCPLS))
      
    }
  }
  
  if (sum(ov_mat_down) > 0){
    pdf("ov_mat_down.pdf")
    NMF::aheatmap(ov_mat_down, Rowv = NA, Colv = NA, color = "Blues", txt = ov_mat_down,
                  main = "Overlap between Giotto findICG and CCPLS")
    dev.off()
  } else {
    print("No ovelapping in down-regulated genes.")
  }
  
}