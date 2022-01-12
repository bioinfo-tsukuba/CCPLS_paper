library(Giotto)

Giotto_preprocess <- function(my_working_dir, coord_file, exp_file, annot_file){
  
  expr_path = fs::path(my_working_dir, exp_file)
  loc_path = fs::path(my_working_dir, coord_file)
  meta_path = fs::path(my_working_dir, annot_file)
  
  ## first merge location and additional metadata
  SS_locations = data.table::fread(loc_path)
  meta_fields = data.table::fread(meta_path)

  SS_expression = data.table::fread(expr_path)
  SS_expression_t <- SS_expression %>% t()
  colnames(SS_expression_t) <- SS_locations$cell_ID
  
  ## create Giotto object
  my_giotto_object <- createGiottoObject(raw_exprs = SS_expression_t,
                                         spatial_locs = SS_locations)
  
  ## add additional annotation if wanted
  my_giotto_object = addCellMetadata(my_giotto_object,
                               new_metadata = meta_fields,
                               by_column = T,
                               column_cell_ID = 'cell_ID')
  
  ## filter
  my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                             expression_threshold = 1,
                             gene_det_in_min_cells = 10,
                             min_det_genes_per_cell = 10,
                             expression_values = c('raw'),
                             verbose = T)
  
  ## normalize
  my_giotto_object <- normalizeGiotto(gobject = my_giotto_object, scalefactor = 6000, verbose = T)
  
  ## add gene & cell statistics
  my_giotto_object <- addStatistics(gobject = my_giotto_object)
  
  ## adjust expression matrix for technical or known variables
  my_giotto_object <- adjustGiottoMatrix(gobject = my_giotto_object, expression_values = c('normalized'),
                                   batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),
                                   return_gobject = TRUE,
                                   update_slot = c('custom'))

  ## highly variable genes (HVG)
  my_giotto_object <- calculateHVG(gobject = my_giotto_object, method = 'cov_loess', difference_in_cov = 0.1)
  
  ## select genes based on HVG and gene statistics, both found in gene metadata
  gene_metadata = fDataDT(my_giotto_object)
  featgenes = gene_metadata[hvg == 'yes' & perc_cells > 4 & mean_expr_det > 0.5]$gene_ID
  
  ## run PCA on expression values (default)
  my_giotto_object <- runPCA(gobject = my_giotto_object, genes_to_use = featgenes, scale_unit = F, center = F)

  ## run UMAP and tSNE on PCA space (default)
  my_giotto_object <- runUMAP(my_giotto_object)

  my_giotto_object <- runtSNE(my_giotto_object)

  ## sNN network (default)
  my_giotto_object <- createNearestNetwork(gobject = my_giotto_object)

  ## cell typeを代入する
  my_giotto_object@cell_metadata$cell_types <- meta_fields[,3]
  
  my_giotto_object = createSpatialNetwork(gobject = my_giotto_object)
  
  browser()
  
  high_expressed_genes = gene_metadata[mean_expr_det > 1.31]$gene_ID
  
  ICGscoresHighGenes =  findICG(gobject = my_giotto_object,
                                selected_genes = high_expressed_genes,
                                spatial_network_name = 'Delaunay_network',
                                cluster_column = 'cell_types',
                                diff_test = 'permutation',
                                adjust_method = 'fdr',
                                nr_permutations = 2000, 
                                do_parallel = T, cores = 4)
  
  plotCellProximityGenes(SS_seqfish, cpgObject = ICGscoresHighGenes, 
                         method = 'dotplot')
  
  cell_proximities = cellProximityEnrichment(gobject = my_giotto_object,
                                             cluster_column = 'cell_types',
                                             spatial_network_name = 'Delaunay_network',
                                             adjust_method = 'fdr',
                                             number_of_simulations = 2000)
  
  cellProximityNetwork(gobject = my_giotto_object,
                       CPscore = cell_proximities,
                       remove_self_edges = T,
                       only_show_enrichment_edges = T)
  
  return(my_giotto_object)

}

#####関数読み込み######
run_findICG <- function(my_giotto_object_for_ICG){
  result <-  findICG(gobject = my_giotto_object_for_ICG,
                     cluster_column = 'cell_types')
  result_2 <- filterICG(result)
  return(result_2)
}



############
# seqFISH
dataset_seqFISH <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/seqFISH+_real_data',
  exp_file_all = c("exp_mat_all.tsv"),
  coord_file_all = c("coord_mat_all.tsv"),
  annot_file_all = c("annot_mat_all.tsv"),
  rdata_path = "seqFISH_findICG_result.Rdata",
  tsv_prefix = "seqFISH_"   
)

# Seq-Scope
dataset_Seq_Scope <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/processed_seq_scope_data',
  exp_file_all = c("exp_mat_all.tsv"),
  coord_file_all = c("coord_mat_all.tsv"),
  annot_file_all = c("annot_mat_all.tsv"),
  rdata_path = "Seq_Scope_findICG_result.Rdata",
  tsv_prefix = "Seq_Scope_"   
)

# toy model 4
dataset_toy_model_4_1 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_1_alpha_1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_1_alpha_1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_1_alpha_1"   
)
dataset_toy_model_4_2 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_1_alpha_0.3_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_1_alpha_0.3_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_1_alpha_0.3"   
)
dataset_toy_model_4_3 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_1_alpha_0.1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_1_alpha_0.1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_1_alpha_0.1"   
)
dataset_toy_model_4_4 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.3_alpha_1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.3_alpha_1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.3_alpha_1"   
)
dataset_toy_model_4_5 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.3_alpha_0.3_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.3_alpha_0.3_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.3_alpha_0.3"   
)
dataset_toy_model_4_6 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.3_alpha_0.1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.3_alpha_0.1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.3_alpha_0.1"   
)
dataset_toy_model_4_7 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.1_alpha_1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.1_alpha_1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.1_alpha_1"   
)
dataset_toy_model_4_8 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.1_alpha_0.3_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.1_alpha_0.3_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.1_alpha_0.3"   
)
dataset_toy_model_4_9 <- list(
  working_dir = '/Users/takaho/Gitclone/CCPLS_Packaging/Work/2021_12_giotto_run/simulated_data_4',
  exp_file_all = c("exp_mat_w_0.1_alpha_0.1_all.tsv"),
  coord_file_all = c("coord_mat_toy.tsv"),
  annot_file_all = c("annot_mat_toy.tsv"),
  rdata_path = "toy_model_4_findICG_w_0.1_alpha_0.1_result.Rdata",
  tsv_prefix = "toy_model_4_matrix_w_0.1_alpha_0.1"   
)