setwd("~/CCPLS_paper/Work/2022_06_giotto/Seq_Scope")
download.file("https://figshare.com/ndownloader/files/36417207",
              "~/CCPLS_paper/Work/2022_06_giotto/Seq_Scope/exp_mat_Seq_Scope.tsv")

library(dplyr)
library(ggplot2)

source("findICG-preprocess-to-run-Seq-Scope.R")
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
  
  res.sel.var <- readRDS("res.sel.var.seq_scope.rds")
  ret_h <- heatmap_ICG(df_ICG, res.sel.var)
  ret_c <- compare_CCPLS(df_ICG, res.sel.var, ret_h)
  
  if (data_flag == "Seq_Scope"){
    
    print(ret_c$ct_cmb[[19]])
    
    bg_go <- unlist(read.table("3_gene_bg_list.txt", header = TRUE))
    # print(ret_c$ct_cmb[[19]])
    icg_19 <- ret_c$ICG_symbol_list_up[[19]]
    
    library(clusterProfiler)
    library(DOSE)
    library(org.Mm.eg.db)
    
    hs <- org.Mm.eg.db
    bg_go_2 <- biomaRt::select(hs,
                               keys = bg_go,
                               columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
    icg_19_2 <- biomaRt::select(hs,
                                keys = icg_19,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
    ego <- clusterProfiler::enrichGO(gene = icg_19_2[,"ENTREZID"],
                                     universe      = bg_go_2[,"ENTREZID"],
                                     OrgDb         = hs,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     readable      = TRUE)
    ego_2 <- simplify(ego)
    if(nrow(ego_2) == 0){
      print("No enriced in B.cell-IgA => B.cell-Immature in ICGs")
    }
  
    ccp_g <- ret_c$CCPLS_symbol_list_up[[19]]

    ## GO results 
    # Tbc1d20/Hnf4a/Sox18/Rapgef2/Nup210l/Myadm/Sidt2/Gpx1/Psap/Pcdh15/Xbp1/Foxa1/Atf4/Fshr/Yipf6
    # go_g <- c("Tbc1d20a", "Hnf4a", "Sox18", "Rapgef2", "Nup210l", "Myadm", "Sidt2",
    #          "Gpx1", "Psap", "Pcdh15", "Xbp1", "Foxa1", "Atf4", "Fshr", "Yipf6")
    # coef_mat_imm <- res.sel.var$sig_coef_mat_list[[3]]
    # go_g_flag <- colnames(coef_mat_imm) %in% go_g
    # coef_go_g <- coef_mat_imm["B.cell-IgA", go_g_flag]
          
    library(readxl)
    ed_g_1 <- read_xlsx("GO_term_summary_20220706_072357.xlsx") # Epithelial cell development
    ed_g_2 <- unique(ed_g_1$Symbol)
  
    library(gplots)
    pdf("venn.pdf")
    data = list(CCPLS = ccp_g, ICG = icg_19, GO = ed_g_2)
    venn(data)
    dev.off()
    
  }
  
  
  # print(str(findICG_result))
  
}

run_findICG_batch(dataset_Seq_Scope, run_opt = TRUE)