setwd("~/CCPLS_paper/Work/2022_06_giotto/seqFISH_plus")

library(dplyr)
library(ggplot2)

source("findICG-preprocess-to-run-seqFISH_plus.R")
source("collectICG.R")
source("heatmap_ICG.R")
source("compare_CCPLS.R")
# source("evaluateICG.R")

run_findICG_batch <- function(dataset, run_opt){
  
  working_dir <- dataset$working_dir
  coord_file_all <- dataset$coord_file_all
  exp_file_all <- dataset$exp_file_all
  annot_file_all <- dataset$annot_file_all
  rdata_path <- dataset$rdata_path
  rdata_path2 <- dataset$rdata_path2
  tsv_prefix <- dataset$tsv_prefix
  data_flag <- dataset$data_flag
  
  print(working_dir)
  print(coord_file_all)
  print(exp_file_all)
  print(annot_file_all)
  print(rdata_path)
  print(rdata_path2)
  print(tsv_prefix)
  
  sample_ids <- gsub(".tsv$", "", gsub("exp_mat_", "", basename(exp_file_all)))
  
  if (run_opt){
    
    findICG_result <- NULL
    findICG_result <- Giotto_preprocess_to_run(working_dir, coord_file_all, exp_file_all, annot_file_all)

    save(findICG_result, file = rdata_path)
    
  } else {
    
    load(rdata_path)
    
  }
  
  df_ICG <- collectICG(tsv_prefix, findICG_result)
  save(df_ICG, file = rdata_path2)

  heatmap_ICG(df_ICG)
  
  res.main <- readRDS("res.main.seqfish_plus.rds")
  compare_CCPLS(df_ICG, res.main)
  
  print(str(findICG_result))
  
}

run_findICG_batch(dataset_seqFISH_plus, run_opt = TRUE)