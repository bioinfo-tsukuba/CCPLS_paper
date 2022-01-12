readSeqfishPlusData <- function(my_working_dir){
  
  dir_exp_mat <- paste(my_working_dir, "/original_data/exp_mat_giotto.txt", sep = "")
  dir_annot_mat <- paste(my_working_dir, "/original_data/annot_mat_giotto.txt", sep = "")
  dir_coord_mat <- paste(my_working_dir, "/original_data/coord_mat_giotto.txt", sep = "")

  exp_mat <- read.table(dir_exp_mat, header = T)
  exp_mat <- as.matrix(exp_mat)
  
  cell_num <- nrow(exp_mat)
  
  cell_ID_column <- matrix(cell_num, nrow = cell_num, ncol = 1)
  
  for (i in 1:cell_num){
    cell_ID_column[i,1] <- paste0("s", i)
  }
  
  rownames(exp_mat) <- cell_ID_column
  
  annot_mat <- read.table(dir_annot_mat, header = T)
  annot_mat[,1] <- as.character(cell_ID_column)
  colnames(annot_mat) <- c("cell_ID", "FOV", "cell_type")
  annot_mat[,3] <- stringr::str_replace_all(annot_mat[,3], " ", ".")
  
  coord_mat <- read.table(dir_coord_mat, header = T)
  coord_mat <- cbind(as.character(cell_ID_column), coord_mat)
  colnames(coord_mat) <- c("cell_ID", "X", "Y")
  
  
  return(list(exp_mat = exp_mat,
              annot_mat = annot_mat,
              corrd_mat = coord_mat))
  
}