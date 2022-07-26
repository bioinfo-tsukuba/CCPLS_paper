Sim_exp_mat_v15.4_internal <- function(exp_mat, coord_mat, annot_mat,
                                    annot_score_mat, gene_score_mat,
                                    HVG_opt, PCA_opt,
                                    HVG_extract_num, PCA_extract_num,
                                    exp_mat_all_sep_list,
                                    fet_mat_sep_list,
                                    coef_sep_list,
                                    e_g_i_const,
                                    weight_const,
                                    my_working_dir, output_dir,
                                    output_opt,
                                    HVG_list_l5,
                                    exp_mat_l5){

  ##############################
  ## exp_mat_sim を生成する関数
  ##############################
  
  ## g: gene g, k: cell type k, i: cell i, j: feature j

  set.seed(123)
  
  library(Seurat)
  library(dplyr)
  library(qvalue)

  ### Get Mean-CV and Do Mean-CV regression
  cell_type_num <- coef_sep_list %>% length()

  weighted_gene_sep_list <- vector("list", cell_type_num)
  
  exp_mat_seur_sep_list <- vector("list", cell_type_num)
  HVG_sep_list <- vector("list", cell_type_num)

  ## u_{g,k} を実データの exp_mat から取得して格納する変数
  mean_list <- vector("list", cell_type_num)

  ## e_{g,i} \sim N(0, sigma) における sigma を実データの exp_mat から取得して格納する変数
  sd_list <- vector("list", cell_type_num)
  CV_list <- vector("list", cell_type_num)
  sd_list_s <- vector("list", cell_type_num)

  ## Mean-CV regression した結果を格納する変数
  fit_list <- vector("list", cell_type_num)
  log_CV_est_list <- vector("list", cell_type_num)
  CV_est_list <- vector("list", cell_type_num)
  
  for (cell_type_ind in 1:1){

    # Convert to Seurat Object
    # seur_obj <- CreateSeuratObject(counts = t(exp_mat_all_sep_list[[cell_type_ind]]))
    seur_obj <- CreateSeuratObject(counts = t(exp_mat_l5))
    
    # Normalize
    # seur_obj <- NormalizeData(seur_obj)
    seur_obj <- NormalizeData(seur_obj)

    # Identify Highly Variable Genes
    # seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = HVG_extract_num)
    seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = HVG_extract_num)
    
    HVG_flag <- seur_obj@assays$RNA@meta.features[,"vst.variable"]
    HVG_list_sim <- rownames(seur_obj@assays$RNA@meta.features)[HVG_flag]
    HVG_sep_list[[cell_type_ind]] <- HVG_list_sim
    sum(HVG_list_sim %in% HVG_list_l5) == 2000

    exp_mat_seur_1 <- seur_obj@assays$RNA@data %>% as.matrix() %>% t()
    exp_mat_seur_2 <- exp_mat_seur_1[, HVG_list_sim]
    exp_mat_seur_sep_list[[cell_type_ind]] <- exp_mat_seur_2
      
    mean_list[[cell_type_ind]] <- apply(exp_mat_seur_sep_list[[cell_type_ind]], 2, mean)

    sd_list[[cell_type_ind]] <- apply(exp_mat_seur_sep_list[[cell_type_ind]], 2, sd)
    
    # CV calculation
    CV_list[[cell_type_ind]] <- apply(exp_mat_seur_sep_list[[cell_type_ind]], 2, sd) / mean_list[[cell_type_ind]]
    
    # Mean-CV regression
    mean_list[[cell_type_ind]] %>% as_tibble() -> mean_df
    colnames(mean_df) <- "Mean"

    CV_list[[cell_type_ind]] %>% as_tibble() -> cv_df
    colnames(cv_df) <- "CV"

    cbind(mean_df, cv_df) %>% as_tibble() -> mean_cv_df

    mean_cv_df %>% mutate(CV2 = CV ^2) %>%
      mutate(log_Mean = log10(Mean)) %>%
      mutate(log_CV = log10(CV)) -> mean_cv_df2

    fit_list[[cell_type_ind]] <- loess(
                                        formula = log_CV ~ log_Mean,
                                        data = mean_cv_df2,
                                        span = 0.3 # default value in Seurat
                                       )

    # Get estimated CV
    log_CV_est_list[[cell_type_ind]] <- fit_list[[cell_type_ind]]$fitted
    CV_est_list[[cell_type_ind]] <- 10 ^ log_CV_est_list[[cell_type_ind]]
    
  }
  
  # 生成する
  ## y'_{g,i} = \sum _j { w_{g,j} x_{i,j} } + u_{g,k} + e_{g,i}

  # Cell X Gene (HVG only)
  exp_mat_raw_sep_list_sim <- vector("list", cell_type_num)
  exp_mat_raw_scaled_sep_list_sim <- vector("list", cell_type_num)
  
  exp_mat_raw_sim_k_list <- list("vector", 1)
  sd_mat_raw_sim_k_list <- list("vector", 1)
  
  for (cell_type_ind in 1:cell_type_num){
    
    HVG_list_k <- HVG_sep_list[[cell_type_ind]]
    
    row_name_k <- exp_mat_all_sep_list[[cell_type_ind]] %>% rownames()
    
    col_name_k <- HVG_sep_list[[cell_type_ind]]
    y_num_k <- HVG_extract_num
    cell_num_k <- exp_mat_all_sep_list[[cell_type_ind]] %>% nrow()
    
    fet_mat_k_raw <- fet_mat_sep_list[[cell_type_ind]]
    neig_ind <- grep("neig", colnames(fet_mat_k_raw))
    fet_mat_k <- scale(fet_mat_k_raw[, neig_ind], center = TRUE, scale = TRUE)
    
    # sd calculation
    sd_k <- CV_est_list[[cell_type_ind]] * mean_list[[cell_type_ind]]
    # CV_k <- CV_est_list[[cell_type_ind]] # 10^ スケール
    CV_k <- log_CV_est_list[[cell_type_ind]]
    names(CV_k) <- names(mean_list[[cell_type_ind]])
    
    w_k_list <- list("vector", 1)
    
    weighted_gene_sep_list <- list("vector", cell_type_num)
    
    non_null_flag <- col_name_k[1] != "NULL" && col_name_k[1] != "No_significant"
    
    if (non_null_flag){
    
      w_k <- coef_sep_list[[cell_type_ind]]
      w_k_list[[1]] <- w_k
      
      weighted_gene_list_k <- c()
      
      for (fet_ind in 1:nrow(w_k)){
        
          weighted_gene_list_k_fet <- names(w_k[fet_ind, (abs(w_k) > 0)[fet_ind, ]])
          add_gene_list <- setdiff(weighted_gene_list_k_fet, weighted_gene_list_k)
        
          weighted_gene_list_k <- append(weighted_gene_list_k, add_gene_list)
        
        }
      
      weighted_gene_sep_list[[cell_type_ind]] <- weighted_gene_list_k
      
    } else {
      
      w_k <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
      
    }
    
    exp_mat_raw_sim_k <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
    colnames(exp_mat_raw_sim_k) <- col_name_k
    rownames(exp_mat_raw_sim_k) <- row_name_k
      
    sd_mat_raw_sim_k <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
    colnames(sd_mat_raw_sim_k) <- col_name_k
    rownames(sd_mat_raw_sim_k) <- row_name_k
      
    w_k <- w_k_list[[1]]
      
    for (cell_ind in 1:cell_num_k){
      
      fet_i <- fet_mat_k[cell_ind,] %>% unlist()
      
      for (y_ind in 1:y_num_k){
        
        target_gene_name <- colnames(w_k)[y_ind]
        
        w_g_k <- w_k[, y_ind] %>% unlist()
        
        e_g_i <- rnorm(1, mean = 0, sd = sd_k[target_gene_name]) * e_g_i_const
        # e_g_i <- rnorm(1, mean = 0, sd = abs(CV_k[target_gene_name])) * e_g_i_const

        if (non_null_flag){
          sum_term <- w_g_k %*% fet_i
        } else {
          sum_term <- 0
        }
        
        # y'_{g,i} = \sum _j { w_{g,j} x_{i,j} } + u_{g,k} + e_{g,i}

        exp_mat_raw_sim_k[cell_ind, y_ind] <- sum_term
        sd_mat_raw_sim_k[cell_ind, y_ind] <- e_g_i
          
      }
    
    }
    
    exp_mat_raw_sim_k_list[[1]] <- exp_mat_raw_sim_k
    sd_mat_raw_sim_k_list[[1]] <- sd_mat_raw_sim_k
    
    apply(exp_mat_raw_sim_k, 2, var)
    
    exp_mat_raw_sim_k_sum <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
    colnames(exp_mat_raw_sim_k_sum) <- col_name_k
    rownames(exp_mat_raw_sim_k_sum) <- row_name_k
    
    exp_mat_raw_sim_k_sum <- exp_mat_raw_sim_k_sum + exp_mat_raw_sim_k_list[[1]] + sd_mat_raw_sim_k_list[[1]]
    
    exp_mat_raw_sep_list_sim[[1]] <- exp_mat_raw_sim_k_sum
    exp_mat_raw_scaled_sep_list_sim[[1]] <- exp_mat_raw_sim_k_sum
    
  }
  
  if (output_opt){
    
    if (sum(exp_mat_raw_sep_list_sim[[1]] != 0)){
      png(paste0(my_working_dir, output_dir, "/exp_mat_scaled.png"))
      NMF::aheatmap(exp_mat_raw_sep_list_sim[[1]])
      dev.off()
      
      png(paste0(my_working_dir, output_dir, "/feature_mat.png"))
      NMF::aheatmap(fet_mat_sep_list[[1]][, grep("neig", colnames(fet_mat_sep_list[[1]]))], Colv = NA)
      dev.off()
      
      png(paste0(my_working_dir, output_dir, "/feature_mat_scaled.png"))
      NMF::aheatmap(scale(fet_mat_sep_list[[1]][, grep("neig", colnames(fet_mat_sep_list[[1]]))]), Colv = NA)
      dev.off()
    } else {
      sink(paste0(my_working_dir, output_dir, "/simulated_exp_mat_scaled.txt"))
      print("sum(exp_mat_raw_sep_list_sim[[1]] is 0.")
      sink()
    }
    
  }
  
  return(list(exp_mat_raw_sep_list_sim = exp_mat_raw_sep_list_sim,
              exp_mat_raw_scaled_sep_list_sim = exp_mat_raw_scaled_sep_list_sim,
              exp_mat_raw_sim_k_list = exp_mat_raw_sim_k_list,
              sd_mat_raw_sim_k_list = sd_mat_raw_sim_k_list,
              fet_mat_sep_list = fet_mat_sep_list,
              mean_list = mean_list,
              sd_list = sd_list,
              CV_list = CV_list,
              CV_est_list = CV_est_list,
              HVG_sep_list = HVG_sep_list,
              weighted_gene_sep_list = weighted_gene_sep_list))

}
