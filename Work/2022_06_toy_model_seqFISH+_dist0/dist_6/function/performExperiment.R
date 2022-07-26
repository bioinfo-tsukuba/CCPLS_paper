performExperiment <- function(my_working_dir,
                              output_dir,
                              weight_const,
                              e_g_i_const,
                              HVG_list_l5,
                              exp_mat_l5){
  
  coef_toy[1,1] <- weight_const
  coef_toy[2,2] <- weight_const
  coef_toy[3,3] <- weight_const
  coef_toy[3,1] <- weight_const
    
  print(coef_toy)
  
  ## Weight preparing
  res.sep.mat.seqFISHplus <- readRDS(paste0(my_working_dir, "/object/res.sep.mat.rds"))
  
  exp_mat_all_sep_list <- vector("list", length = 1)
  HVG_list <- HVG_list_l5
  exp_mat_all_sep_list[[1]] <- res.sep.mat.seqFISHplus$exp_mat_orig_sep_list[[10]][,HVG_list] # L5.eNeuron
  
  # print(res.sep.mat.seqFISHplus$cell_type_list[[10]])
  
  source(paste0(my_working_dir, "/function/prepareWeight2.R"))
  res.prepareWeight <- prepareWeight2(my_working_dir, coef_toy, HVG_list, weight_const)
  
  coef_sep_list <- res.prepareWeight$coef_sep_list
  weighted_SRG_sep_list <- res.prepareWeight$weighted_SRG_sep_list
  
  # Convert annotation of original cell types to A, B, D and D
  source(paste0(my_working_dir, "/function/convertAnnotMatv13.1.R"))
  annot_mat_v13.1 <- convertAnnotMatV13.1(annot_mat_seqFISHplus)
  cell_type_list_v13.1_pre <- unique(annot_mat_v13.1[, "cell_type"])
  cell_type_list_v13.1 <- cell_type_list_v13.1_pre[order(cell_type_list_v13.1_pre)]
  
  annot_score_mat_v13.1 <- matrix(0, nrow = nrow(annot_mat_v13.1),
                                  ncol = length(cell_type_list_v13.1))
  rownames(annot_score_mat_v13.1) <- rownames(annot_mat_v13.1)
  colnames(annot_score_mat_v13.1) <- cell_type_list_v13.1
  
  for (i in 1:length(cell_type_list_v13.1)){
    each_cell_type <- cell_type_list_v13.1[i]
    annot_score_mat_v13.1[annot_mat_v13.1[, "cell_type"] == each_cell_type, i] <- 1
  }
  
  # Make feature matrix of toy model v13.1
  fet_mat_sep_list <- vector("list", length = 1)
  fet_mat_sep_list_2 <- vector("list", length = 1)
  
  A_flag <- annot_mat_v13.1[, "cell_type"] == "A"

  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegMakeFet.R"))
  res.make.features <- cellCellRegMakeFet(coord_mat = coord_mat_seqFISHplus,
                                          annot_mat = annot_mat_v13.1)
  
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegMakeFet_2.R"))
  res.make.features.2 <- cellCellRegMakeFet_2(coord_mat = coord_mat_seqFISHplus,
                                              annot_mat = annot_mat_v13.1,
                                              dist_diff = 6)
  
  
  fet_mat_sep_list[[1]] <- res.make.features$fet_mat_norm[A_flag,]
  fet_mat_sep_list_2[[1]] <- res.make.features.2$fet_mat_norm[A_flag,]
  
  png("fet_mat_2.png")
  NMF::aheatmap(res.make.features.2$fet_mat_norm[A_flag, 5:8],
                Rowv = NA, Colv = NA)
  dev.off()
  
  ## Prepare arguments
  exp_mat <- exp_mat_seqFISHplus
  coord_mat <- coord_mat_seqFISHplus
  annot_mat <- annot_mat_v13.1
  annot_score_mat <- res.make.features$annot_score_mat
  gene_score_mat <- gene_score_mat_seqFISHplus
  
  # experimental condition
  HVG_opt = TRUE
  PCA_opt = FALSE
  HVG_extract_num = 2000
  PCA_extract_num = 10
  perm_num = 100
  species_in_your_data = c("human", "mouse")[2]
  dev_opt = "pls"
  
  # Simulate
  source(paste0(my_working_dir, "/function/Sim_exp_mat_v15.4_internal.R"))
  res.sim <- Sim_exp_mat_v15.4_internal(exp_mat, coord_mat, annot_mat,
                                        annot_score_mat, gene_score_mat,
                                        HVG_opt, PCA_opt,
                                        HVG_extract_num, PCA_extract_num,
                                        exp_mat_all_sep_list,
                                        fet_mat_sep_list,
                                        coef_sep_list, e_g_i_const,
                                        weight_const,
                                        my_working_dir, output_dir,
                                        output_opt = TRUE,
                                        HVG_list_l5,
                                        exp_mat_l5)
  
  exp_mat_raw_sep_list_sim <- res.sim$exp_mat_raw_sep_list_sim
  fet_mat_sep_list <- res.sim$fet_mat_sep_list
  cell_type_list <- "A"
  
  res.sep.mat <- list(exp_mat_sep_list = res.sim$exp_mat_raw_sep_list_sim,
                      # fet_mat_sep_list = res.sim$fet_mat_sep_list)
                      fet_mat_sep_list = fet_mat_sep_list_2)
  
  #### Perform CCPLS
  
  ### Estimate by PLS regression
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegEst.R"))
  res.estimate <- cellCellRegEst(fet_mat_orig = fet_mat_sep_list[[1]],
                                 coord_mat = coord_mat,
                                 annot_score_mat = annot_score_mat,
                                 exp_mat_sep_list = res.sep.mat$exp_mat_sep_list,
                                 fet_mat_sep_list = res.sep.mat$fet_mat_sep_list,
                                 estimate_cell_type_opt = "all",
                                 cell_type_list = cell_type_list,
                                 dev_opt = dev_opt)
  
  res.estimate$cell_type_list <- "A"
  
  set.seed(1)
  
  ### Select significant variable
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegSelVar.R"))
  res.sel.var <- cellCellRegSelVar(res.estimate, res.sep.mat, HVG_extract_num, dev_opt = "kmeans")
  
  ### Get heatmap object
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegHeatmap.R"))
  res.heatmap <- cellCellRegHeatmap(res.sel.var, output_dir)
  
  ### Get bipartite graph object
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegGraph.R"))
  res.graph <- cellCellRegGraph(res.sel.var, output_dir)
  
  # ARI, sensitivity and specificity
  ans_clust <- rep("non-SRG", HVG_extract_num)
  names(ans_clust) <- colnames(res.sep.mat$exp_mat_sep_list[[1]])
  ans_clust[names(ans_clust) %in% weighted_SRG_sep_list[[1]]] <- "SRG1"
  ans_clust[names(ans_clust) %in% weighted_SRG_sep_list[[2]]] <- "SRG2"
  ans_clust[names(ans_clust) %in% weighted_SRG_sep_list[[3]]] <- "SRG3"
  ans_clust[names(ans_clust) %in% weighted_SRG_sep_list[[4]]] <- "non-SRG"
  
  bin_mat <- res.sel.var$sig_coef_mat_bin_2_list[[1]]
  
  if (!is.null(bin_mat)){
    
    if (bin_mat != "NULL"){
      
      SRG1 <- c()
      SRG2 <- c()
      SRG3 <- c()
      non_SRG <- c()
      
      for (gene_id in 1:ncol(bin_mat)){
        
        t_gene <- colnames(bin_mat)[gene_id]
        t_bin <- bin_mat[,gene_id]
        SRG1_judge <- t_bin["neig_A"] == 1 & t_bin["neig_B"] == 0 & t_bin["neig_C"] == 1 & t_bin["neig_D"] == 0
        SRG2_judge <- t_bin["neig_A"] == 0 & t_bin["neig_B"] == 1 & t_bin["neig_C"] == 0 & t_bin["neig_D"] == 0
        SRG3_judge <- t_bin["neig_A"] == 0 & t_bin["neig_B"] == 0 & t_bin["neig_C"] == 1 & t_bin["neig_D"] == 0
        non_SRG_judge <- t_bin["neig_A"] == 0 & t_bin["neig_B"] == 0 & t_bin["neig_C"] == 0 & t_bin["neig_D"] == 0
        if(SRG1_judge){SRG1 <- append(SRG1, t_gene)}
        if(SRG2_judge){SRG2 <- append(SRG2, t_gene)}
        if(SRG3_judge){SRG3 <- append(SRG3, t_gene)}
        if(non_SRG_judge){non_SRG <- append(non_SRG, t_gene)}
        
      }
      
      est_clust <- rep("Others", HVG_extract_num)
      names(est_clust) <- colnames(res.sep.mat$exp_mat_sep_list[[1]])
      
      cluster_vec <- res.sel.var$gene_cluster_vec_list[[1]]
      cluster_num <- max(cluster_vec)
      
      for (cluster_ind in 1:cluster_num){
        cluster_gene_name <- names(cluster_vec[cluster_vec == cluster_ind])
        cluster_gene_num <- length(cluster_gene_name)
        
        SRG1_flag <- (sum(SRG1 %in% cluster_gene_name) / cluster_gene_num) > 0.5
        SRG2_flag <- (sum(SRG2 %in% cluster_gene_name) / cluster_gene_num) > 0.5
        SRG3_flag <- (sum(SRG3 %in% cluster_gene_name) / cluster_gene_num) > 0.5
        non_SRG_flag <- (sum(non_SRG %in% cluster_gene_name) / cluster_gene_num) > 0.5
        
        if (SRG1_flag){
          est_clust[cluster_gene_name] <- "SRG1"
        } else if (SRG2_flag){
          est_clust[cluster_gene_name] <- "SRG2"
        } else if (SRG3_flag){
          est_clust[cluster_gene_name] <- "SRG3"
        } else if (non_SRG_flag){
          est_clust[cluster_gene_name] <- "non_SRG"
        }
        
        sink(paste0(my_working_dir, output_dir, "/classify_cluster_", cluster_ind, "_.txt"))
        print(paste0("SRG1 flag is ", SRG1_flag))
        print(paste0("SRG2 flag is ", SRG2_flag))
        print(paste0("SRG3 flag is ", SRG3_flag))
        print(paste0("non-SRG flag is ", non_SRG_flag))
        sink()
        
      }
      
      non_SRG_cluster <- setdiff(colnames(exp_mat_all_sep_list[[1]]), names(res.sel.var$gene_cluster_vec_list[[1]]))
      non_SRG_cluster_gene_num <- length(non_SRG_cluster)
      
      SRG1_flag <- (sum(SRG1 %in% non_SRG_cluster) / non_SRG_cluster_gene_num) > 0.5
      SRG2_flag <- (sum(SRG2 %in% non_SRG_cluster) / non_SRG_cluster_gene_num) > 0.5
      SRG3_flag <- (sum(SRG3 %in% non_SRG_cluster) / non_SRG_cluster_gene_num) > 0.5
      non_SRG_flag <- (sum(non_SRG %in% non_SRG_cluster) / non_SRG_cluster_gene_num) > 0.5

      if (SRG1_flag){
        est_clust[non_SRG_cluster] <- "SRG1"
      } else if (SRG2_flag){
        est_clust[non_SRG_cluster] <- "SRG2"
      } else if (SRG3_flag){
        est_clust[non_SRG_cluster] <- "SRG3"
      } else if (non_SRG_flag){
        est_clust[non_SRG_cluster] <- "non_SRG"
      }
      
      ari <- mclust::adjustedRandIndex(ans_clust, est_clust)
      
      sink(paste0(my_working_dir, output_dir, "/ARI.txt"))
      print(paste0("ARI is ", ari))
      sink()
      
      sig_coef_mat <- res.sel.var$sig_coef_mat_with_zero_list[[1]]
      ans_coef_mat <- coef_sep_list[[1]]
      pcc_all <- cor.test(sig_coef_mat, ans_coef_mat)$estimate
      pcc_p <- cor.test(sig_coef_mat, ans_coef_mat)$p.value
      
      Others <- names(est_clust[est_clust == "Others"])
      
      ########## Correspondence table #############################3
      #                     Actual
      #                   SRG1, SRG2, SRG3, non-SRG
      # Estimated   SRG1  c11   c12   c13   c14    
      #             SRG2  c21   c22   c23   c24
      #             SRG3  c31   c32   c33   c34
      #          non-SRG  c41   c42   c43   c44
      #           Others  c51   c52   c53   c54
      #
      
      c11 <- sum(SRG1 %in% weighted_SRG_sep_list[[1]])
      c12 <- sum(SRG1 %in% weighted_SRG_sep_list[[2]])
      c13 <- sum(SRG1 %in% weighted_SRG_sep_list[[3]])
      c14 <- sum(SRG1 %in% weighted_SRG_sep_list[[4]])
      
      c21 <- sum(SRG2 %in% weighted_SRG_sep_list[[1]])
      c22 <- sum(SRG2 %in% weighted_SRG_sep_list[[2]])
      c23 <- sum(SRG2 %in% weighted_SRG_sep_list[[3]])
      c24 <- sum(SRG2 %in% weighted_SRG_sep_list[[4]])
      
      c31 <- sum(SRG3 %in% weighted_SRG_sep_list[[1]])
      c32 <- sum(SRG3 %in% weighted_SRG_sep_list[[2]])
      c33 <- sum(SRG3 %in% weighted_SRG_sep_list[[3]])
      c34 <- sum(SRG3 %in% weighted_SRG_sep_list[[4]])
      
      c41 <- sum(non_SRG %in% weighted_SRG_sep_list[[1]])
      c42 <- sum(non_SRG %in% weighted_SRG_sep_list[[2]])
      c43 <- sum(non_SRG %in% weighted_SRG_sep_list[[3]])
      c44 <- sum(non_SRG %in% weighted_SRG_sep_list[[4]])
      
      c51 <- sum(Others %in% weighted_SRG_sep_list[[1]])
      c52 <- sum(Others %in% weighted_SRG_sep_list[[2]])
      c53 <- sum(Others %in% weighted_SRG_sep_list[[3]])
      c54 <- sum(Others %in% weighted_SRG_sep_list[[4]])
      
      c_tab <- rbind(c(c11, c12, c13, c14),
                     c(c21, c22, c23, c24),
                     c(c31, c32, c33, c34),
                     c(c41, c42, c43, c44),
                     c(c51, c52, c53, c54))
      
      tp_srg1 <- c11
      fp_srg1 <- c12 + c13 + c14
      fn_srg1 <- c21 + c31 + c41 + c51
      tn_srg1 <- sum(c_tab) - tp_srg1 - fp_srg1 - fn_srg1
      
      prec_srg1 <- tp_srg1 / (tp_srg1 + fp_srg1)
      rec_srg1 <- tp_srg1 / (tp_srg1 + fn_srg1)
      
      tp_srg2 <- c22
      fp_srg2 <- c21 + c23 + c24
      fn_srg2 <- c12 + c32 + c42 + c52
      tn_srg2 <- sum(c_tab) - tp_srg2 - fp_srg2 - fn_srg2
      
      prec_srg2 <- tp_srg2 / (tp_srg2 + fp_srg2)
      rec_srg2 <- tp_srg2 / (tp_srg2 + fn_srg2)
      
      tp_srg3 <- c33
      fp_srg3 <- c31 + c32 + c34
      fn_srg3 <- c13 + c23 + c43 + c53
      tn_srg3 <- sum(c_tab) - tp_srg3 - fp_srg3 - fn_srg3
      
      prec_srg3 <- tp_srg3 / (tp_srg3 + fp_srg3)
      rec_srg3 <- tp_srg3 / (tp_srg3 + fn_srg3)
      
      
      tp_non_srg <- c44
      fp_non_srg <- c41 + c42 + c43
      fn_non_srg <- c14 + c24 + c34 +c54
      tn_non_srg <- sum(c_tab) - tp_non_srg - fp_non_srg - fn_non_srg
      
      prec_non_srg <- tp_non_srg / (tp_non_srg + fp_non_srg)
      rec_non_srg <- tp_non_srg / (tp_non_srg + fn_non_srg)
      
      sink(paste0(my_working_dir, output_dir, "/precision_recall.txt"))
      print(paste0("PCC is ", pcc_all))
      print("In SRG1,")
      print(paste0("Precision is ", prec_srg1))
      print(paste0("Recall is ", rec_srg1))
      print("")
      print("In SRG2,")
      print(paste0("Precision is ", prec_srg2))
      print(paste0("Recall is ", rec_srg2))
      print("")
      print("In SRG3,")
      print(paste0("Precision is ", prec_srg3))
      print(paste0("Recall is ", rec_srg3))
      print("")
      print("In non-SRG,")
      print(paste0("Precision is ", prec_non_srg))
      print(paste0("Recall is ", rec_non_srg))
      sink()
      
      ind_df <- as_tibble(rbind(ari, pcc_all,
                                prec_srg1, prec_srg2, prec_srg3, prec_non_srg,
                                rec_srg1, rec_srg2, rec_srg3, rec_non_srg))
      colnames(ind_df) <- "Value"
      
      ind_list <- c("Adjusted rand index", "Pearson correlation coefficient",
                    "Precision of cluster 1", "Precision of cluster 2", "Precision of cluster 3", "Precision of cluster 4",
                    "Recall of cluster 1", "Recall of cluster 2", "Recall of cluster 3", "Recall of cluster 4")
      
      ind_df %>% mutate(Index = ind_list) -> ind_df_2
      
      p_ind <- ggplot(ind_df_2, aes(x = Index, y = Value)) +
        geom_point(size = 3) +
        ylim(0, 1) +
        theme_bw() +
        xlab("") +
        theme(text = element_text(size = 24)) +
        theme(#panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line(colour = "grey60", linetype = "dashed"),
              axis.text.x = element_text(angle = 60, hjust = 1))
      
      ggsave(paste0(my_working_dir, output_dir, "/index.png"),
             plot = p_ind, width = 8)

      ggsave(paste0(my_working_dir, output_dir, "/index.pdf"),
             plot = p_ind, width = 8)
      
      ind_df_2[1:2, ] %>% mutate(Pseudo = "Pseudo") -> ind_df_3 # "Adjusted rand index" and "Pearson correlation coefficient"
      
      p_ind_2 <- ggplot(ind_df_3, aes(x = Index, y = Value, fill = Pseudo)) +
        geom_col() +
        ylim(0, 1) +
        theme_bw() +
        scale_fill_brewer(palette = "Set1") +
        xlab("") +
        theme(text = element_text(size = 12)) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.position = 'none')
      
      ggsave(paste0(my_working_dir, output_dir, "/index_2.png"),
             plot = p_ind_2, width = 2)
      
      ggsave(paste0(my_working_dir, output_dir, "/index_2.pdf"),
             plot = p_ind_2, width = 2)
      
      
    } else {
      
      sink(paste0(my_working_dir, output_dir, "/index.txt"))
      print("No feature and gene are detected.")
      sink()
      
      pcc_all <- 0
      pcc_p <- NULL
      ari <- 0
      prec_srg1 <- 0
      prec_srg2 <- 0
      prec_srg3 <- 0
      prec_non_srg <- 0
      rec_srg1 <- 0
      rec_srg2 <- 0
      rec_srg3 <- 0
      rec_non_srg <- 0
      
    }
    
  } else {
    
    pcc_all <- 0
    pcc_p <- NULL
    ari <- 0
    prec_srg1 <- 0
    prec_srg2 <- 0
    prec_srg3 <- 0
    prec_non_srg <- 0
    rec_srg1 <- 0
    rec_srg2 <- 0
    rec_srg3 <- 0
    rec_non_srg <- 0
    
  }

  
  var_tot <- apply(res.sim$exp_mat_raw_sep_list[[1]], 2, var)
  var_sd <- apply(res.sim$sd_mat_raw_sim_k_list[[1]], 2, var)
  var_prop <- (var_tot - var_sd) / var_tot
  var_prop_mean <- mean(var_prop)
  
  sink(paste0(my_working_dir, output_dir, "/variance_proportion.txt"))
  print(paste0("Variance proportion is ", var_prop_mean))
  sink()
  
  return(list(ari = ari,
              pcc_all = pcc_all,
              pcc_p = pcc_p,
              prec_srg1 = prec_srg1,
              prec_srg2 = prec_srg2,
              prec_srg3 = prec_srg3,
              prec_non_srg = prec_non_srg,
              rec_srg1 = rec_srg1,
              rec_srg2 = rec_srg2,
              rec_srg3 = rec_srg3,
              rec_non_srg = rec_non_srg,
              var_prop_mean = var_prop_mean))
  
}