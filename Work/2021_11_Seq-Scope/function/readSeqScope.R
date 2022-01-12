readSeqScope <- function(file_list, file_index){
  
  library("Seurat")
  
  output_dir <- paste0("/results_v0.1.8_", file_list[file_index])
  data_obj <- readRDS(paste0(my_working_dir, "/original_data/", file_list[file_index]))
  
  cell_num <- nrow(data_obj@meta.data)
  gene_num <- nrow(data_obj@assays$Spatial@counts)
  
  exp_mat <- t(as.matrix(data_obj@assays$Spatial@counts))
  
  coord_mat <- matrix(0, nrow = cell_num, ncol = 3)
  colnames(coord_mat) <- c("cell_ID", "x_coord", "y_coord")
  coord_mat[,"x_coord"] <- data_obj$X
  coord_mat[,"y_coord"] <- data_obj$Y
  coord_mat[,"cell_ID"] <- colnames(data_obj@assays$Spatial@counts)
  
  annot_mat <- matrix(0, nrow = cell_num, ncol = 3)
  colnames(annot_mat) <- c("cell_ID", "FOV", "cell_types")
  
  annot_mat[,"cell_types"] <- as.character(data_obj@active.ident)
  annot_mat[,"cell_types"] <- stringr::str_replace_all(annot_mat[,"cell_types"], " ", ".")
  annot_mat[,"cell_ID"] <- colnames(data_obj@assays$Spatial@counts)
  
  
  return(list(exp_mat = exp_mat,
              coord_mat = coord_mat,
              annot_mat = annot_mat))
  
}