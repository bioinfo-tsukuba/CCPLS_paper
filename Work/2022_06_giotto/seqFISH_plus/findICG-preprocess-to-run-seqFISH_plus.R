library(Giotto)

Giotto_preprocess_to_run <- function(my_working_dir, coord_file, exp_file, annot_file){
  
  createGiottoInstructions()
  
  expr_path = fs::path(my_working_dir, exp_file)
  loc_path = fs::path(my_working_dir, coord_file)
  meta_path = fs::path(my_working_dir, annot_file)
  
  ## first merge location and additional metadata
  SS_locations <- data.table::fread(loc_path)
  meta_fields <- data.table::fread(meta_path)

  SS_expression <- data.table::fread(expr_path)
  SS_expression_t <- SS_expression %>% t()
  colnames(SS_expression_t) <- SS_locations$cell_ID
  
  ## create Giotto object
  my_giotto_object <- createGiottoObject(raw_exprs = SS_expression_t,
                                         spatial_locs = SS_locations)
  
  ## add additional annotation if wanted
  my_giotto_object <- addCellMetadata(my_giotto_object,
                               new_metadata = meta_fields,
                               by_column = T,
                               column_cell_ID = 'cell_ID')
  
  ## normalize
  my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
  
  ## add gene & cell statistics
  my_giotto_object <- addStatistics(gobject = my_giotto_object)
  
  ## adjust expression matrix for technical or known variables
  my_giotto_object <- adjustGiottoMatrix(gobject = my_giotto_object)

  ## cell typeを代入する
  my_giotto_object@cell_metadata$cell_types <- meta_fields[,3]
  
  my_giotto_object <- createSpatialNetwork(gobject = my_giotto_object)
  
  ICGscoresHighGenes <- findICG(gobject = my_giotto_object, cluster_column = 'cell_types')
  
  result <- filterICG(ICGscoresHighGenes)
  
  return(result)

}


############
# seqFISH plus
dataset_seqFISH_plus <- list(
  working_dir = "~/CCPLS_paper/Work/2022_06_giotto/seqFISH_plus",
  exp_file_all = c("exp_mat_seqFISH_plus.tsv"),
  coord_file_all = c("coord_mat_seqFISH_plus.tsv"),
  annot_file_all = c("annot_mat_seqFISH_plus.tsv"),
  rdata_path = "seqFISH_plus_result.Rdata",
  rdata_path2 = "seqFISH_plus_result_2.Rdata",
  tsv_prefix = "seqFISH_plus",
  data_flag = "seqFISH_plus"
)