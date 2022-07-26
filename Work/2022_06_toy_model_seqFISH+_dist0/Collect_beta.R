Collect_beta <- function(beta_vec, beta_dir){
  
  exp_v <- c()
  ari_v <- c()
  pcc_v <- c()
  prec_srg1_v <- c()
  prec_srg2_v <- c()
  prec_srg3_v <- c()
  prec_non_srg_v <- c()
  rec_srg1_v <- c()
  rec_srg2_v <- c()
  rec_srg3_v <- c()
  rec_non_srg_v <- c()
  
  for (i in 1:length(beta_vec)){
    
    res.list <- readRDS(paste0("~/CCPLS_Paper/Work/2022_06_toy_model_seqFISH+_dist0/dist_",
                                beta_dir[i], "/results_collect/res.list.rds"))
    res.list.1 <- res.list[[9]] # w = 1, alpha = 1
    
    exp_v <- append(exp_v, paste0(beta_vec[i]))
    
    ari_v <- append(ari_v, res.list.1$ari)
    pcc_v <- append(pcc_v, res.list.1$pcc_all)
    prec_srg1_v <- append(prec_srg1_v, res.list.1$prec_srg1)
    prec_srg2_v <- append(prec_srg2_v, res.list.1$prec_srg2)
    prec_srg3_v <- append(prec_srg3_v, res.list.1$prec_srg3)
    prec_non_srg_v <- append(prec_non_srg_v, res.list.1$prec_non_srg)
    rec_srg1_v <- append(rec_srg1_v, res.list.1$prec_srg1)
    rec_srg2_v <- append(rec_srg2_v, res.list.1$prec_srg2)
    rec_srg3_v <- append(rec_srg3_v, res.list.1$prec_srg3)
    rec_non_srg_v <- append(rec_non_srg_v, res.list.1$prec_non_srg)
    
  }
  
  # Beta
  rep(exp_v, 10, each = 1) -> exp_v_2
  
  # Index
  rep("Adjusted rand index", length(beta_vec)) %>%
    append(rep("Pearson correlation coefficient", length(beta_vec))) %>%
    append(rep("Precision of cluster 1", length(beta_vec))) %>%
    append(rep("Precision of cluster 2", length(beta_vec))) %>%
    append(rep("Precision of cluster 3", length(beta_vec))) %>%
    append(rep("Precision of cluster 4", length(beta_vec))) %>%
    append(rep("Recall of cluster 1", length(beta_vec))) %>%
    append(rep("Recall of cluster 2", length(beta_vec))) %>%
    append(rep("Recall of cluster 3", length(beta_vec))) %>%
    append(rep("Recall of cluster 4", length(beta_vec))) -> ind_v
    
  # Value
  ari_v %>% append(pcc_v) %>%
    append(prec_srg1_v) %>%
    append(prec_srg2_v) %>%
    append(prec_srg3_v) %>%
    append(prec_non_srg_v) %>%
    append(rec_srg1_v) %>%
    append(rec_srg2_v) %>%
    append(rec_srg3_v) %>%
    append(rec_non_srg_v) -> val_v
  val_v[is.na(val_v)] <- 0
  
  # Combine
  as_tibble(exp_v_2) %>% cbind(ind_v) %>% cbind(val_v) -> df_1
  colnames(df_1) <- c("Beta", "Index", "Value")  
  group_by(df_1, "Index")
  df_1[,"Beta"] <- as.double(df_1[,"Beta"])
  
  p <- ggplot(df_1, aes(x = Beta, y = Value, color = Index)) +
          geom_line() +
          scale_x_log10() +
          scale_color_brewer(palette = "Paired") +
          theme_bw() +
          ggtitle("Performance of ten indexes against dist0 perturbation") +
          facet_wrap(. ~ Index, ncol = 3)
  ggsave("perturb.png", p, width = 10, height = 8)
  ggsave("perturb.pdf", p, width = 10, height = 8)
  
}

setwd("~/CCPLS_Paper/Work/2022_06_toy_model_seqFISH+_dist0")
Collect_beta(beta_vec = c(10, 6, 3, 1, 0.75, 0.5, 1/3, 0.2),
             beta_dir = c("10", "6", "3", "1", "0.75", "0.5", "0.33", "0.2"))