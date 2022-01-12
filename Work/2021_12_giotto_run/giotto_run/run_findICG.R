setwd("~/CCPLS_paper/Work/2021_12_giotto_run/giotto_run")

library(dplyr)
library(ggplot2)

source("findICG-preprocess-to-run.R")
source("collectICG.R")
source("evaluateICG.R")

run_findICG_batch <- function(dataset, run_opt){
  
  working_dir <- dataset$working_dir
  coord_file_all <- dataset$coord_file_all
  exp_file_all <- dataset$exp_file_all
  annot_file_all <- dataset$annot_file_all
  rdata_path <- dataset$rdata_path
  tsv_prefix <- dataset$tsv_prefix
  data_flag <- dataset$data_flag
  
  print(working_dir)
  print(coord_file_all)
  print(exp_file_all)
  print(annot_file_all)
  print(rdata_path)
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
  
  if (data_flag == "toy_w_1_alpha_1"){
    res_eval <- evaluate_ICG_toy(df_ICG)
    compare_res_1(res_eval)

  } else if (data_flag == "Seq-Scope") {
    print(evaluate_ICG_seqscope(df_ICG))
  }
  
  print(str(findICG_result))
  
}

run_findICG_batch(dataset_toy_model_4_1, run_opt = TRUE)