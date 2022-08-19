Sim_exp_mat_rev_internal <- function(exp_mat, coord_mat, annot_mat,
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
  
  ## exp_mat_l5とHVG_list_l5からノイズ行列を生成する
  var_v <- apply(exp_mat_l5[, HVG_list_l5], 2, var)
  mean_v <- apply(exp_mat_l5[, HVG_list_l5], 2, mean)
  
  # shape = v / u^2
  shape_v <- var_v / (mean_v)^2
  # scale = u / v
  scale_v <- mean_v / var_v
  
  noise_mat_raw <- matrix(0, nrow = nrow(exp_mat_l5), ncol = 2000)
  rownames(noise_mat_raw) <- rownames(exp_mat_l5)
  colnames(noise_mat_raw) <- HVG_list_l5
  
  for (h_i in 1:2000){
    noise_h_i <- rgamma(nrow(exp_mat_l5),
                        shape = shape_v[h_i],
                        scale = scale_v[h_i])
    noise_mat_raw[, h_i] <- noise_h_i
  }
  
  ## Apply noise_mat_raw to Seurat
  
  library(Seurat)
  library(dplyr)
  library(qvalue)

  # Apply Seurat to noise_mat_raw
  seur_obj_sim <- CreateSeuratObject(counts = t(noise_mat_raw))
  seur_obj_sim <- NormalizeData(seur_obj_sim)
  seur_obj_sim@assays$RNA@data %>% as.matrix() %>% t() -> noise_mat
  noise_mat_z <- scale(noise_mat)
  
  # Apply Seurat to exp_mat_l5
  #seur_obj <- CreateSeuratObject(counts = t(exp_mat_l5))
  #seur_obj <- NormalizeData(seur_obj)
  #seur_obj@assays$RNA@data %>% as.matrix() %>% t() -> exp_mat_s
  
  # sim 
  cell_type_num <- coef_sep_list %>% length()
  weighted_gene_sep_list <- vector("list", cell_type_num)
  
  exp_mat_seur_sep_list <- vector("list", cell_type_num)
  HVG_sep_list <- vector("list", cell_type_num)
  
  # 生成する
  ## y'_{g,i} = \sum _j { w_{g,j} x_{i,j} } + u_{g,k} + e_{g,i}

  # Cell X Gene (HVG only)
  exp_mat_raw_sep_list_sim <- vector("list", cell_type_num)
  exp_mat_raw_scaled_sep_list_sim <- vector("list", cell_type_num)
  
  exp_mat_raw_sim_k_list <- list("vector", 1)
  noise_mat_raw_sim_k_list <- list("vector", 1)
  
  for (ct_i in 1:cell_type_num){
    
    HVG_list_k <- HVG_list_l5
    row_name_k <- rownames(exp_mat_l5)
    col_name_k <- HVG_list_l5
    y_num_k <- 2000
    cell_num_k <- nrow(exp_mat_l5)
    
    fet_mat_k_raw <- fet_mat_sep_list[[ct_i]]
    neig_ind <- grep("neig", colnames(fet_mat_k_raw))
    fet_mat_k <- scale(fet_mat_k_raw[, neig_ind], center = TRUE, scale = TRUE)
    
    w_k_list <- list("vector", 1)
    
    weighted_gene_sep_list <- list("vector", cell_type_num)
    
    non_null_flag <- col_name_k[1] != "NULL" && col_name_k[1] != "No_significant"
    
    if (non_null_flag){
    
      w_k <- coef_sep_list[[ct_i]]
      w_k_list[[1]] <- w_k
      
      weighted_gene_list_k <- c()
      
      for (fet_ind in 1:nrow(w_k)){
        
          weighted_gene_list_k_fet <- names(w_k[fet_ind, (abs(w_k) > 0)[fet_ind, ]])
          add_gene_list <- setdiff(weighted_gene_list_k_fet, weighted_gene_list_k)
        
          weighted_gene_list_k <- append(weighted_gene_list_k, add_gene_list)
        
        }
      
      weighted_gene_sep_list[[ct_i]] <- weighted_gene_list_k
      
    } else {
      
      w_k <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
      
    }
    
    exp_mat_raw_sim_k <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
    colnames(exp_mat_raw_sim_k) <- col_name_k
    rownames(exp_mat_raw_sim_k) <- row_name_k
      
    w_k <- w_k_list[[1]]
      
    for (cell_ind in 1:cell_num_k){
      
      fet_i <- fet_mat_k[cell_ind,] %>% unlist()
      
      for (y_ind in 1:y_num_k){
        
        w_g_k <- w_k[, y_ind] %>% unlist()
        
        if (non_null_flag){
          sum_term <- w_g_k %*% fet_i
        } else {
          sum_term <- 0
        }
        
        # y'_{g,i} = \sum _j { w_{g,j} x_{i,j} } + u_{g,k} + e_{g,i}

        exp_mat_raw_sim_k[cell_ind, y_ind] <- sum_term

      }
    
    }
    
    exp_mat_raw_sim_k_list[[1]] <- exp_mat_raw_sim_k
    noise_mat_raw_sim_k_list[[1]] <- noise_mat_z
    
    exp_mat_raw_sim_k_sum <- matrix(0, nrow = cell_num_k, ncol = y_num_k)
    colnames(exp_mat_raw_sim_k_sum) <- col_name_k
    rownames(exp_mat_raw_sim_k_sum) <- row_name_k
    
    exp_mat_raw_sim_k_sum <- exp_mat_raw_sim_k_sum + exp_mat_raw_sim_k_list[[1]] + (noise_mat_raw_sim_k_list[[1]] * e_g_i_const)
    
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
              sd_mat_raw_sim_k_list = noise_mat_raw_sim_k_list,
              fet_mat_sep_list = fet_mat_sep_list,
              # mean_list = mean_list,
              # sd_list = sd_list,
              # CV_list = CV_list,
              # CV_est_list = CV_est_list,
              HVG_sep_list = HVG_sep_list,
              weighted_gene_sep_list = weighted_gene_sep_list))

}
