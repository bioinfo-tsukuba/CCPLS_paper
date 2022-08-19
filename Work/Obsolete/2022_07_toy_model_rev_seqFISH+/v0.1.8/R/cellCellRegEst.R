#' Perform PLS regression
#' 
#' @param data_4_pls Data object for PLS regression
#' @param cv_opt Option of cross-validation type
#' @param cv_num Option of cross-validation number
#' @param component_num Number of component
#' 
#' @importFrom pls plsr
#' 
buildModel <- function(data_4_pls,
                       cv_opt = c("CV", "LOOCV")[1],
                       cv_num = 10,
                       component_num = NULL){
  
  if (cv_opt == "CV"){
    
    # cv_num fold CV
    
    if (is.null(component_num)){
      res_pls <- plsr(gene ~ feature, data = data_4_pls,
                      ncomp = qr(data_4_pls$feature)$rank,
                      scale = TRUE,
                      validation = "CV", segments = cv_num)
    } else {
      res_pls <- plsr(gene ~ feature, data = data_4_pls,
                      ncomp = component_num,
                      scale = TRUE,
                      validation = "CV", segments = cv_num)
    }
    
    RMSE_cv_all <- sqrt(res_pls$validation$PRESS / nrow(data_4_pls$feature))
    RMSE_cv_all_sum <- colSums(RMSE_cv_all)
    
    opt_comp <- RMSE_cv_all_sum[RMSE_cv_all_sum == min(RMSE_cv_all_sum)] %>% names()
    opt_comp_num <- strsplit(opt_comp, split = " ")[[1]][1] %>% as.double()
    
  } else if (cv_opt == "LOOCV") {
    
    # LOOCV
    
    if (is.null(component_num)){
      res_pls <- plsr(gene ~ feature, data = data_4_pls,
                      ncomp = qr(data_4_pls$feature)$rank,
                      scale = TRUE,
                      validation = "LOO")
    } else {
      res_pls <- plsr(gene ~ feature, data = data_4_pls,
                      ncomp = componetn_num,
                      scale = TRUE,
                      validation = "LOO")
    }
    
    RMSE_cv_all <- sqrt(res_pls$validation$PRESS / nrow(data_4_pls$feature))
    RMSE_cv_all_sum <- colSums(RMSE_cv_all)
    
    opt_comp <- RMSE_cv_all_sum[RMSE_cv_all_sum == min(RMSE_cv_all_sum)] %>% names()
    opt_comp_num <- strsplit(opt_comp, split = " ")[[1]][1] %>% as.double()
    
  }
  
  return(list(res_pls = res_pls,
              opt_comp_num = opt_comp_num,
              RMSE_cv_all = RMSE_cv_all,
              RMSE_cv_all_sum = RMSE_cv_all_sum))
  
}

#' Estimate PLS regression model
#'
#' @param fet_mat_orig Feature matrix
#' @param coord_mat Coordinate matrix
#' @param annot_score_mat Annotation score matrix
#' @param exp_mat_sep_list Separated expression matrix as list
#' @param fet_mat_sep_list Separated feature matrix as list
#' @param estimate_cell_type_opt Option of selecting target cell types in estimation
#' @param cell_type_list Cell type list
#' @param perm_opt Option of permutation
#' @param component_num_list List of component number
#' @param dev_opt Option for development mode
#'
#' @export
#' 
cellCellRegEst <- function(fet_mat_orig, coord_mat,
                           annot_score_mat,
                           exp_mat_sep_list,
                           fet_mat_sep_list,
                           estimate_cell_type_opt = "all", # estimate_cell_type_opt = cell_type_model_SG
                           cell_type_list,
                           perm_opt = FALSE,
                           component_num_list = NULL,
                           dev_opt = c("ropls", "pls")[1]
){
  
  ### For development
  library(tictoc) # TO BE DELETED with tic() toc()
  
  start_time <- Sys.time()
  
  ## PLS library
  library(ropls)
  library(pls)
  library(parallel)
  
  # prepare function
  beera <- function(expr){
    tryCatch(expr,
             error = function(e){
               # message("An error occurred:\n", e)
               return("An error occurred")
             })
  }
  
  target_cell_type_index <- 0
  gene_list <- list()
  
  ## choose column index for estimating model
  if (estimate_cell_type_opt == "all"){
    
    cell_type_model_SG <- c()
    
    # "self-"
    self_col_index <- grep("self", colnames(fet_mat_orig))
    estimate_col_index <- seq(1:length(cell_type_list))
    
    component_num_list <- rep(0, length = length(estimate_col_index))
    estimate_cell_type_list <- cell_type_list
    
  } else {
    
    cell_type_model_SG <- estimate_cell_type_opt
    
    # "self-"
    self_col_index <- grep("self", colnames(fet_mat_orig))
    estimate_col_index <- rep(0, length(cell_type_model_SG))
    
    for(estimate_num in 1:length(cell_type_model_SG)){
      estimate_col_index[estimate_num] <- grep(cell_type_model_SG[estimate_num], cell_type_list)
    }
    
    estimate_cell_type_list <- cell_type_list[estimate_col_index]
    
  }
  
  estimate_col_num <- length(estimate_col_index)
  self_cell_type_list <- cell_type_list
  
  CCPLS_result <- vector("list", estimate_col_num)
  
  # "self-" ごとに PLS regression
  for (i in 1:estimate_col_num){
    
    print(paste0("Under calculation in ", self_cell_type_list[estimate_col_index[i]], "..."))

    # "self-"
    fet_mat_sep <- as.matrix(fet_mat_sep_list[[estimate_col_index[i]]])
    exp_mat_sep <- as.matrix(exp_mat_sep_list[[estimate_col_index[i]]])
    
    fet_mat_sep_rem <- fet_mat_sep[, - self_col_index]
    
    # Remove zero variance
    fet_mat_sep_rem <- fet_mat_sep_rem[,apply(fet_mat_sep_rem, 2, var) != 0]
    
    fet_mat_sep_rem_norm <- scale(fet_mat_sep_rem, center = TRUE, scale = TRUE) %>% as.data.frame()
    exp_mat_sep_norm <- scale(exp_mat_sep, center = TRUE, scale = TRUE) %>% as.data.frame()
    
    # plsr
    tic()
        
    data_4_pls <- list()
    data_4_pls$feature <- as.matrix(fet_mat_sep_rem_norm)
    data_4_pls$gene<- as.matrix(exp_mat_sep_norm)
        
    if ((length(data_4_pls$feature) != 0) && (!is.na(data_4_pls$gene)[1]) && (colnames(data_4_pls$gene) != "NULL")){
        
      returned_value <- buildModel(data_4_pls = data_4_pls,
                                   cv_opt = c("CV", "LOOCV")[1],
                                   cv_num = 10)
          
      component_num_list[i] <- returned_value$opt_comp_num
          
    } else {
          
      returned_value <- NULL
          
      CCPLS_result[[i]] <- list(paste0("No model was built."),
                                paste0("Model content is empty."))
          
    }
        
    toc()
      
    if ((length(data_4_pls$feature != 0)) && (!is.na(data_4_pls$gene)[1]) && (colnames(data_4_pls$gene) != "NULL")){
      target_cell_type_index <- target_cell_type_index + 1
        
      CCPLS_result[[i]] <- list(paste0(self_cell_type_list[estimate_col_index[i]]),
                                returned_value)
        
      gene_list[[target_cell_type_index]] <- colnames(exp_mat_sep_norm)
        
    } else {
        
      CCPLS_result[[i]] <- list(paste0("No model was built."),
                                paste0("Model content is empty."))
        
    }
      
  }
    
  for (j in 1:length(estimate_cell_type_list)){
    
    judge_model_estimated <- base::grep(estimate_cell_type_list[j], CCPLS_result[[j]][[1]]) # [[j]][[1]] に細胞型の名前が入っている
    
    if (length(judge_model_estimated) != 0){
      cell_type_model_SG <- base::append(cell_type_model_SG, self_cell_type_list[j])
    }
    
  }
  
  finish_time <- Sys.time()
  
  
  return(list(CCPLS_result = CCPLS_result,
              self_cell_type_list = self_cell_type_list,
              cell_type_model_SG = cell_type_model_SG,
              cell_type_list = cell_type_list,
              gene_list = gene_list,
              component_num_list = component_num_list,
              start_time = start_time,
              finish_time = finish_time))
  
}