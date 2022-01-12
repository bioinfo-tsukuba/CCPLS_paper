evaluate_ICG_toy <- function(df_ICG){
 
  answer_dir <- c("answer_gene")
  
  cluster_list <- c("cluster_1.txt",
                    "cluster_2.txt",
                    "cluster_3.txt",
                    "cluster_4.txt")

  all_gene <- c()
  for (cluster_id in 1:length(cluster_list)){
    all_gene <- append(all_gene,
                       unlist(read.table(paste0(answer_dir, "/", cluster_list[cluster_id]))))
  }
   
  df_ICG_A <- df_ICG[df_ICG[, "receiver_cell_type"] == "A", ]
  df_ICG_A_HVG <- df_ICG_A[df_ICG_A$sig_gene_symbol %in% all_gene, ]
  df_ICG_A_HVG_upreg <- df_ICG_A_HVG[df_ICG_A_HVG$log2fc > 0, ]
  df_ICG_A_HVG_downreg <- df_ICG_A_HVG[df_ICG_A_HVG$log2fc < 0, ]
  
  cluster_1_flag <- df_ICG_A_HVG_upreg$sender_cell_type == "A" | df_ICG_A_HVG_upreg$sender_cell_type == "C"
  cluster_2_flag <- df_ICG_A_HVG_upreg[, "sender_cell_type"] == "B"
  cluster_3_flag <- df_ICG_A_HVG_upreg[, "sender_cell_type"] == "C"
  others_flag <- df_ICG_A_HVG_upreg[, "sender_cell_type"] == "D"
  
  est_list <- vector("list", length(cluster_list))
  
  est_list[[1]] <- unlist(df_ICG_A_HVG_upreg[cluster_1_flag, "sig_gene_symbol"])
  est_list[[2]] <- unlist(df_ICG_A_HVG_upreg[cluster_2_flag, "sig_gene_symbol"])
  est_list[[3]] <- unlist(df_ICG_A_HVG_upreg[cluster_3_flag, "sig_gene_symbol"])
  est_list[[4]] <- setdiff(all_gene, df_ICG_A_HVG[,"sig_gene_symbol"])
  est_list[[5]] <- unlist(rbind(df_ICG_A_HVG_upreg[others_flag, "sig_gene_symbol"],
                                df_ICG_A_HVG_downreg[, "sig_gene_symbol"]))
  
  ans_list <- vector("list", length(cluster_list))
  
  ans_list[[1]] <- unlist(read.table(paste0(answer_dir, "/", cluster_list[1])))
  ans_list[[2]] <- unlist(read.table(paste0(answer_dir, "/", cluster_list[2])))
  ans_list[[3]] <- unlist(read.table(paste0(answer_dir, "/", cluster_list[3])))
  ans_list[[4]] <- unlist(read.table(paste0(answer_dir, "/", cluster_list[4])))
  
  ########## Correspondence table #############################3
  #                     Actual
  #                    1     2     3     4
  # Estimated      1  c11   c12   c13   c14    
  #                2  c21   c22   c23   c24
  #                3  c31   c32   c33   c34
  #                4  c41   c42   c43   c44
  #                5  c51   c52   c53   c54
  #
  
  c11 <- sum(est_list[[1]] %in% ans_list[[1]])
  c12 <- sum(est_list[[1]] %in% ans_list[[2]])
  c13 <- sum(est_list[[1]] %in% ans_list[[3]])
  c14 <- sum(est_list[[1]] %in% ans_list[[4]])

  c21 <- sum(est_list[[2]] %in% ans_list[[1]])
  c22 <- sum(est_list[[2]] %in% ans_list[[2]])
  c23 <- sum(est_list[[2]] %in% ans_list[[3]])
  c24 <- sum(est_list[[2]] %in% ans_list[[4]])

  c31 <- sum(est_list[[3]] %in% ans_list[[1]])
  c32 <- sum(est_list[[3]] %in% ans_list[[2]])
  c33 <- sum(est_list[[3]] %in% ans_list[[3]])
  c34 <- sum(est_list[[3]] %in% ans_list[[4]])

  c41 <- sum(est_list[[4]] %in% ans_list[[1]])
  c42 <- sum(est_list[[4]] %in% ans_list[[2]])
  c43 <- sum(est_list[[4]] %in% ans_list[[3]])
  c44 <- sum(est_list[[4]] %in% ans_list[[4]])
  
  c51 <- sum(est_list[[5]] %in% ans_list[[1]])
  c52 <- sum(est_list[[5]] %in% ans_list[[2]])
  c53 <- sum(est_list[[5]] %in% ans_list[[3]])
  c54 <- sum(est_list[[5]] %in% ans_list[[4]])
  
  c_tab <- rbind(c(c11, c12, c13, c14),
                 c(c21, c22, c23, c24),
                 c(c31, c32, c33, c34),
                 c(c41, c42, c43, c44),
                 c(c51, c52, c53, c54))
  
  tp_cluster_1 <- c11
  fp_cluster_1 <- c12 + c13 + c14
  fn_cluster_1 <- c21 + c31 + c41 + c51
  tn_cluster_1 <- sum(c_tab) - tp_cluster_1 - fp_cluster_1 - fn_cluster_1
  
  prec_cluster_1 <- tp_cluster_1 / (tp_cluster_1 + fp_cluster_1)
  rec_cluster_1 <- tp_cluster_1 / (tp_cluster_1 + fn_cluster_1)
  
  tp_cluster_2 <- c22
  fp_cluster_2 <- c21 + c23 + c24
  fn_cluster_2 <- c12 + c32 + c42 + c52
  tn_cluster_2 <- sum(c_tab) - tp_cluster_2 - fp_cluster_2 - fn_cluster_2
  
  prec_cluster_2 <- tp_cluster_2 / (tp_cluster_2 + fp_cluster_2)
  rec_cluster_2 <- tp_cluster_2 / (tp_cluster_2 + fn_cluster_2)
  
  tp_cluster_3 <- c33
  fp_cluster_3 <- c31 + c32 + c34
  fn_cluster_3 <- c13 + c23 + c43 + c53
  tn_cluster_3 <- sum(c_tab) - tp_cluster_3 - fp_cluster_3 - fn_cluster_3
  
  prec_cluster_3 <- tp_cluster_3 / (tp_cluster_3 + fp_cluster_3)
  rec_cluster_3 <- tp_cluster_3 / (tp_cluster_3 + fn_cluster_3)
  
  
  tp_cluster_4 <- c44
  fp_cluster_4 <- c41 + c42 + c43
  fn_cluster_4 <- c14 + c24 + c34 + c54
  tn_cluster_4 <- sum(c_tab) - tp_cluster_4 - fp_cluster_4 - fn_cluster_4
  
  prec_cluster_4 <- tp_cluster_4 / (tp_cluster_4 + fp_cluster_4)
  rec_cluster_4 <- tp_cluster_4 / (tp_cluster_4 + fn_cluster_4)
  
  return(list(prec_cluster_1 = prec_cluster_1,
              prec_cluster_2 = prec_cluster_2,
              prec_cluster_3 = prec_cluster_3,
              prec_cluster_4 = prec_cluster_4,
              rec_cluster_1 = rec_cluster_1,
              rec_cluster_2 = rec_cluster_2,
              rec_cluster_3 = rec_cluster_3,
              rec_cluster_4 = rec_cluster_4))
   
}

compare_res_1 <- function(res_eval){
  
  # Read CCPLS result
  load("~/CCPLS_paper/Work/2021_12_giotto_run/giotto_run/CCPLS_res/res.list.RData")
  
  df_col <- append(as_tibble(res_eval),
                   as_tibble(res.list[[9]])[4:11])
  
  df_col_2 <- as_tibble(unlist(df_col))
  names(df_col_2) <- "Value"
  
  index_list <- append(paste0("Preicision of cluster ", seq(1:4)),
                       paste0("Recall of cluster ", seq(1:4)))
  index_list_2 <- append(index_list, index_list)
  
  method_list <- append(rep("Giotto findICG", 8),
                        rep("CCPLS", 8))
  
  df_col_2 %>% mutate(Index = index_list_2) %>%
    mutate(Method = method_list) -> df_col_3
  
  p <- ggplot(df_col_3, aes(x = Index, y = Value, fill = Method)) +
          geom_col(position = "dodge") +
          theme_bw() +
          scale_fill_brewer(palette = "Set1") +
          theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1)) +
          theme(axis.text.y = element_text(size = 15)) +
          theme(axis.title.y = element_text(size = 15)) +
          xlab("") +
          ggtitle("Comparison with Giotto findICG")

  ggsave("~/CCPLS_paper/Work/2021_12_giotto_run/giotto_run/Figure/comparison_with_Giotto_findICG.pdf",
         width = 10, plot = p)
  
  return()
    
}

evaluate_ICG_seqscope <- function(df_ICG){
  
  df_ICG_b_immat <- df_ICG[df_ICG$receiver_cell_type == "B.cell-Immature", ]
  df_ICG_b_immat_b_iga <- df_ICG_b_immat[df_ICG_b_immat$sender_cell_type == "B.cell-IgA", ]
  df_ICG_b_immat_b_iga_upreg <- df_ICG_b_immat_b_iga[df_ICG_b_immat_b_iga$log2fc > 0, ]
  
  target_gene <- df_ICG_b_immat_b_iga_upreg$sig_gene_symbol
  bg <- unlist(read.table("back_ground_gene_Seq-Scope_B_cell_immature/3_gene_bg_list.txt",
               header = TRUE))
  
  go_opt <- "BP"
  cutoff <- 0.05
  
  ## GO
  library(clusterProfiler)
  library(DOSE)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
  hs <- org.Mm.eg.db
  
  target_gene_2 <- biomaRt::select(hs,
                                   keys = target_gene,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  bg_2 <- biomaRt::select(hs,
                          keys = bg,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")
  
  ego <- clusterProfiler::enrichGO(gene = target_gene_2[,"ENTREZID"],
                                   universe      = bg_2[,"ENTREZID"],
                                   OrgDb         = hs,
                                   ont           = go_opt,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = cutoff,
                                   readable      = TRUE)

  return(paste0("nrow(ego) is ", nrow(ego)))
  
}