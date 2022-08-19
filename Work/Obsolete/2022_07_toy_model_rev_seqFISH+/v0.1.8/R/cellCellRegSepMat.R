#' Perform Seurat
#' 
#' @param exp_mat Expression matrix
#' @param HVG_extract_num Number of HVG extraction
#' @param PCA_opt Option of PCA
#' @param PCA_extract_num Number of PCA extraction
#' 
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat JackStraw
#' @importFrom Seurat ScoreJackStraw
#' @importFrom Seurat Loadings
#' 
#' @importFrom qvalue qvalue
#' 
performSeurat <- function(exp_mat,
                          HVG_extract_num){
  
  library(Seurat)
  library(qvalue)
  
  ### Convert to Seurat Object
  seur_obj <- CreateSeuratObject(counts = t(exp_mat))
  
  ### Normalize
  seur_obj <- NormalizeData(seur_obj)
  
  ### Identify Highly Variable Genes
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = HVG_extract_num)
  
  HVG_flag <- seur_obj@assays$RNA@meta.features[,"vst.variable"]
  HVG_list <- rownames(seur_obj@assays$RNA@meta.features)[HVG_flag]
  
  ### Scaling the data
  all.genes <- rownames(t(exp_mat))
  seur_obj <- ScaleData(seur_obj)
  
  ### Prepare expression matrix of HVG
  exp_mat_seur <- matrix(0, nrow = nrow(exp_mat), ncol = HVG_extract_num)
  exp_mat_seur <- seur_obj@assays$RNA@data[HVG_list,] %>% as.matrix() %>% t()

  return(list(exp_mat_seur = exp_mat_seur,
              HVG_list = HVG_list))
  
}

#' Perform separation of input matrices
#'
#' @param exp_mat Expression matrix
#' @param fet_mat Feature matrix
#' @param annot_score_mat Annotation score matrix
#' @param gene_score_mat Gene score matrix
#' @param HVG_opt Option of HVG extraction
#' @param PCA_opt Option of PCA
#' @param HVG_extract_num Number of HVG extraction
#' @param PCA_extract_num Number of PCA extraction
#' @param my_working_dir Working directory
#'
#' @export
#'
cellCellRegSepMat <- function(exp_mat = exp_mat,
                              fet_mat = fet_mat,
                              annot_score_mat = annot_score_mat,
                              HVG_opt = HVG_opt,
                              HVG_extract_num = HVG_extract_num){
  
  
  print(paste0("=== CCPLS_sep_mat Started...!! ", Sys.time(), " ==="))
  start_time <- Sys.time()
  
  exp_mat <- as.matrix(exp_mat)
  
  cell_row_full_list <- rownames(exp_mat)
  rownames(fet_mat) <- cell_row_full_list
  
  cell_type_list <- colnames(annot_score_mat)
  
  row_list <- rownames(exp_mat)
  col_list <- colnames(exp_mat)
  
  col_list_fet <- colnames(fet_mat)
  
  exp_mat_sep_list <- vector("list", length = ncol(annot_score_mat))
  exp_mat_sep_each <- matrix(0, nrow = nrow(exp_mat), ncol = ncol(exp_mat))
  fet_mat_sep_list <- vector("list", length = ncol(annot_score_mat))
  exp_mat_orig_sep_list <- vector("list", length = ncol(annot_score_mat))
  HVG_sep_list <- vector("list", length = ncol(annot_score_mat))
  
  # FIXME：replace for loop by matrix operation or apply family
  for (cell_type_index in 1:ncol(annot_score_mat)){
    
    exp_mat_sep_each <- exp_mat[annot_mat[,3] == cell_type_list[cell_type_index],]

    # Remove zero expression
    exp_mat_sep_each_del <- exp_mat_sep_each[apply(exp_mat_sep_each, 1, sum) != 0,]
    
    # Remove zero variance
    exp_mat_sep_each_del2 <- exp_mat_sep_each_del[,apply(exp_mat_sep_each_del, 2, var) != 0]
    
    cell_row_list_del2 <- rownames(exp_mat_sep_each_del2)
    fet_mat_sep_each_unnorm <- fet_mat[cell_row_list_del2,]
    fet_mat_sep_each <- scale(fet_mat_sep_each_unnorm, center = TRUE, scale = TRUE)
    
    null_flag <- is.null(nrow(exp_mat_sep_each_del2)) || nrow(exp_mat_sep_each_del2) == 0
    
    if (null_flag){
      
      exp_mat_sep_each_del2 <- matrix(0, nrow = 2, ncol = ncol(exp_mat))
      rownames(exp_mat_sep_each_del2) <- paste0("s", seq(1,2))
      colnames(exp_mat_sep_each_del2) <- col_list
      
      fet_mat_sep_each <- matrix(0, nrow = 2, ncol = ncol(fet_mat))
      colnames(fet_mat_sep_each) <- col_list_fet
      rownames(fet_mat_sep_each) <- paste0("s", seq(1,2))
      
      exp_mat_sep_each_orig <- exp_mat_sep_each_del2
      
    }
    
    if (HVG_opt == TRUE){
      
      if (null_flag){
        
        null_row_names <- rownames(exp_mat_sep_each_del2)
        exp_mat_sep_each_del2 <- matrix(0, nrow = nrow(exp_mat_sep_each_del2), ncol = 1)
        rownames(exp_mat_sep_each_del2) <- null_row_names
        colnames(exp_mat_sep_each_del2) <- "NULL"
        
        loading_mat <- matrix(0, nrow = nrow(exp_mat_sep_each_del2), ncol = 1)
        rownames(loading_mat) <- null_row_names
        colnames(loading_mat) <- "NULL"
        
        exp_mat_sep_each_orig <- exp_mat_sep_each_del2
        
        f.loading_mat <- "NULL"
        q.val.f.loading_mat <- "NULL"
        
        HVG_list <- "NULL"
        
      } else {
        
        exp_mat_sep_each_orig <- exp_mat_sep_each_del
        
        res.perform.Seurat <- performSeurat(exp_mat = exp_mat_sep_each_del2,
                                            HVG_extract_num = HVG_extract_num)
        
        exp_mat_sep_each_del2 <- res.perform.Seurat$exp_mat_seur
        HVG_list <- res.perform.Seurat$HVG_list
        
      }
      
    }
    
    fet_mat_sep_list[[cell_type_index]] <- fet_mat_sep_each
    exp_mat_sep_list[[cell_type_index]] <- exp_mat_sep_each_del2
    exp_mat_orig_sep_list[[cell_type_index]] <- exp_mat_sep_each_orig
    HVG_sep_list[[cell_type_index]] <- HVG_list
    
  }
  
  finish_time <- Sys.time()
  print(paste0("=== Finished...!! ", Sys.time(), " ==="))
  print(paste0("=========="))
  
  
  return(list(exp_mat_sep_list = exp_mat_sep_list,
              fet_mat_sep_list = fet_mat_sep_list,
              exp_mat_orig_sep_list = exp_mat_orig_sep_list,
              HVG_sep_list = HVG_sep_list,
              cell_type_list = cell_type_list,
              start_time = start_time,
              finish_time = finish_time))
  
}

performSeurat <- function(exp_mat,
                          HVG_extract_num){
  
  library(Seurat)
  library(qvalue)
  
  ### Convert to Seurat Object
  seur_obj <- CreateSeuratObject(counts = t(exp_mat))
  
  ### Normalize
  seur_obj <- NormalizeData(seur_obj)
  
  ### Identify Highly Variable Genes
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = HVG_extract_num)
  
  HVG_flag <- seur_obj@assays$RNA@meta.features[,"vst.variable"]
  HVG_list <- rownames(seur_obj@assays$RNA@meta.features)[HVG_flag]
  
  ### Scaling the data
  all.genes <- rownames(t(exp_mat))
  seur_obj <- ScaleData(seur_obj)
  
  ### Prepare expression matrix with HVG
  exp_mat_seur <- matrix(0, nrow = nrow(exp_mat), ncol = HVG_extract_num)
  exp_mat_seur <- seur_obj@assays$RNA@data[HVG_list,] %>% as.matrix() %>% t()
    
  return(list(exp_mat_seur = exp_mat_seur,
              HVG_list = HVG_list))
  
}