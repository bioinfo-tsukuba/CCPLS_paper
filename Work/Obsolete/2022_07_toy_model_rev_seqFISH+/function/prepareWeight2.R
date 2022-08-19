prepareWeight2 <- function(my_working_dir, coef_toy, HVG_list,
                           weight_const){
  
  library(dplyr)
  
  HVG_num <- length(HVG_list)
  
  fet_name_list <- rownames(coef_toy)
  fet_num <- length(fet_name_list)
  
  weighted_SRG_sep_list <- vector("list", length = 4)
  coef_sep_list <- vector("list", length = 1)
  
  ## 重みのデザインと準備
  for (cell_type_ind in 1:1){
    
    coef_each <- matrix(0, nrow = fet_num, ncol = HVG_num)
    
    rownames(coef_each) <- fet_name_list
    colnames(coef_each) <- HVG_list
    
    coef_each["A", 1:500] <- weight_const
    coef_each["B", 501:1000] <- weight_const
    coef_each["C", 1001:1500] <- weight_const

    coef_each["C", 1:500] <- weight_const
    
    coef_sep_list[[cell_type_ind]] <- coef_each
    
    weighted_SRG_sep_list[[1]] <- names(coef_each["A", 1:500]) # SRG1
    weighted_SRG_sep_list[[2]] <- names(coef_each["B", 501:1000]) # SRG2
    weighted_SRG_sep_list[[3]] <- names(coef_each["C", 1001:1500]) # SRG3
    weighted_SRG_sep_list[[4]] <- names(coef_each["D", 1501:2000]) # non_SRG
    
    # Save gene name
    for (gene_id in 1:length(weighted_SRG_sep_list)){
      write.table(weighted_SRG_sep_list[[gene_id]],
                  paste0("~/CCPLS_paper/Work/2021_11_toy_model_seqFISH+/answer_gene/cluster_", gene_id, ".txt"),
                  row.names = FALSE, col.names = FALSE)
    }
    
  }
    

  return(list(coef_sep_list = coef_sep_list,
              weighted_SRG_sep_list = weighted_SRG_sep_list))
  
}