#' Calculate feature matrix
#'
#' @param coord_mat Coordinate matrix
#' @param annot_score_mat Annotation score matrix
#' @param fet_type_opt Option for specifying feature matrix type
#'
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' 
cellCellRegMakeFet_2 <- function(coord_mat, annot_mat, dist_diff){
  
  start_time <- Sys.time()
  print(paste0("=== CCPLS_make_features Started...!! ", Sys.time(), " ==="))
  
  library("stringr")
  
  cell_num <- nrow(coord_mat)
  
  dist_mat <- as.matrix(stats::dist(coord_mat[,c(2,3)]))
  colnames(dist_mat) <- seq(1:ncol(dist_mat))
  rownames(dist_mat) <- seq(1:nrow(dist_mat))
  diag(dist_mat) <- 1000000
  
  ident_cell_type <- unique(annot_mat[,3])[order(unique(annot_mat[,3]))]
  
  annot_score_mat <- matrix(0, nrow = cell_num, ncol = length(ident_cell_type))
  colnames(annot_score_mat) <- ident_cell_type
  cell_type_num <- length(ident_cell_type)
  
  for (cell_id in 1:cell_num){
    ct_judge <- annot_mat[cell_id, 3] == ident_cell_type
    annot_score_mat[cell_id, ct_judge] <- 1
  }

  fet_mat <- matrix(0, cell_num, cell_type_num * 2)
  
  colnames_list <- append(paste0("self_", ident_cell_type),
                          paste0("neig_", ident_cell_type))
  colnames_list <- str_replace_all(colnames_list, " ", "_")
  
  colnames(fet_mat) <- colnames_list
  
  dist0 <- dist_mat %>% min()
  dist0 <- dist0 * dist_diff # For sensitivity analysis
  print(paste0("dist0 is... ", dist0))
  
  neig_col_start <- cell_type_num + 1
  neig_col_end <- cell_type_num * 2
  
  annot_score_mat <- as.matrix(annot_score_mat)
  
  fet_mat[, 1:cell_type_num] <- annot_score_mat
  
  dist_weight_mat <- exp( - dist_mat / dist0)
  
  fet_mat[, neig_col_start:neig_col_end] <-
    fet_mat[, neig_col_start:neig_col_end] + dist_weight_mat %*% annot_score_mat
  
  # Remove zero feature
  fet_mat <- as.data.frame(fet_mat)
  fet_mat <- fet_mat[,colSums(fet_mat)!=0]
  
  fet_mat_norm <- scale(fet_mat, center = T, scale = T)
  
  finish_time <- Sys.time()
  print(paste0("=== Finished...!! ", Sys.time(), " ==="))
  print(paste0("=========="))
  
  return(list(fet_mat_norm = fet_mat_norm,
              fet_mat_unnorm = fet_mat,
              annot_score_mat = annot_score_mat,
              start_time = start_time,
              finish_time = finish_time))
  
}