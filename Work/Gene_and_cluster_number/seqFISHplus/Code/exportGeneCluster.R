exportGeneCluster <- function(res.sel.var, output_dir){
  
  ct_l <- res.sel.var$cell_type_list
  ct_l_2 <- gsub("/", "", ct_l)
  ct_num <- length(ct_l_2)
  
  for (i in 1:ct_num){
    
    gc_vec <- res.sel.var$gene_cluster_vec_list[[i]]
    
    if (gc_vec[1] != "NULL"){
      
      fn <- paste0(output_dir, "Receiver_cell_type_", ct_l_2[i], ".csv")
      
      gc_vec_2 <- gc_vec[order(gc_vec)]
      gc_table <- cbind(names(gc_vec_2), gc_vec_2)
      colnames(gc_table) <- c("Gene symbol", "Cluster number")
      
      write.csv(gc_table, fn, row.names = FALSE)
      
    }
    
  }
  
}

# Run
res.sel.var <- readRDS("~/CCPLS_paper/Work/Gene_and_cluster_number/seqFISHplus/Code/res.sel.var.rds")
output_dir <- "~/CCPLS_paper/Work/Gene_and_cluster_number/seqFISHplus/"
exportGeneCluster(res.sel.var, output_dir)