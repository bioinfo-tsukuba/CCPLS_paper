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
  
  res.sel.var <- readRDS("res.sel.var.seqfish_plus.rds")
  ret_h <- heatmap_ICG(df_ICG, res.sel.var)
  ret_c <- compare_CCPLS(df_ICG, res.sel.var, ret_h)
  
  if (data_flag == "seqFISH_plus"){
    
    # 134/143/144
    print(ret_c$ct_cmb[[134]])
    print(ret_c$ct_cmb[[143]])
    print(ret_c$ct_cmb[[144]])
    
    bg_go <- unlist(read.table("12_gene_bg_list.txt", header = TRUE))

    ast_opc_g <- ret_c$ICG_symbol_list_up[[134]]
    olig_opc_g <- ret_c$ICG_symbol_list_up[[143]]
    opc_opc_g <- ret_c$ICG_symbol_list_up[[144]]
    icg_go <- c(ast_opc_g, olig_opc_g, opc_opc_g)
    
    library(clusterProfiler)
    library(DOSE)
    library(org.Mm.eg.db)
    
    hs <- org.Mm.eg.db
    bg_go_2 <- biomaRt::select(hs,
                               keys = bg_go,
                               columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
    icg_go_2 <- biomaRt::select(hs,
                                keys = icg_go,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
    ego <- clusterProfiler::enrichGO(gene = icg_go_2[,"ENTREZID"],
                                     universe      = bg_go_2[,"ENTREZID"],
                                     OrgDb         = hs,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     readable      = TRUE)
    ego_2 <- simplify(ego)
    if(nrow(ego_2) == 0){
      print("No enriced in astrocyte, Olig and OPC => OPC in ICGs")
    }
    
    # Abl1/Tmem98/Mag
    ccp_g_1 <- ret_c$CCPLS_symbol_list_up[[134]]
    ccp_g_2 <- ret_c$CCPLS_symbol_list_up[[143]]
    ccp_g_3 <- ret_c$CCPLS_symbol_list_up[[144]]
    ccp_g <- union(union(ccp_g_1, ccp_g_2), ccp_g_3)
    # coef_mat_opc <- res.sel.var$sig_coef_mat_list[[12]]
    # go_g_flag <- colnames(coef_mat_opc) %in% go_g
    # coef_go_g <- coef_mat_opc[, go_g_flag]
    
    library(readxl)
    gd_g_1 <- read_xlsx("GO_term_summary_20220706_060244.xlsx") # glial cell differentiations
    gd_g_2 <- unique(gd_g_1$Symbol)
    
    library(gplots)
    pdf("venn.pdf")
    data = list(CCPLS = ccp_g, ICG = icg_go, GO = gd_g_2)
    venn(data)
    dev.off()
    
  }
  
  # print(str(findICG_result))
  
}

run_findICG_batch(dataset_seqFISH_plus, run_opt = TRUE)