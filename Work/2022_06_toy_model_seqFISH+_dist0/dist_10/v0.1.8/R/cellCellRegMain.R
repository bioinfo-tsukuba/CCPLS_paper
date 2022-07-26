#' Perform CCPLS
#' 
#' Arguments for cellCellRegressionMakeFeature
#' @param coord_mat Coordinate matrix
#' @param annot_mat Annotaiton matrix
#' 
#' Arguments for cellCellRegressionEstimate
#' @param exp_mat_sep_list Separated expression matrix as list
#' @param fet_mat_sep_list Separated feature matrix as list
#' @param cell_type_list Cell type list
#' @param dev_opt Option for development mode
#' @param HVG_opt Option for extract HVG
#' @param PCA_opt Option for PCA execution
#' @param HVG_extract_num Specification of HVG extraction number
#' @param PCA_extract_num Specification of PCA extraction number
#' @param sim_opt Option for simulated data usage
#' 
#' Arguments for permutation test
#' @param perm_num Trial number of permutation test
#' 
#' @export
#' 
cellCellRegMain <- function(exp_mat = exp_mat,
                            coord_mat = coord_mat,
                            annot_mat = annot_mat,
                            HVG_opt = TRUE,
                            HVG_extract_num = 2000,
                            dev_opt = "pls",
                            my_working_dir,
                            output_dir){
  

  ### Make feature matrix
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegMakeFet.R"))
  res.make.features <- cellCellRegMakeFet(coord_mat = coord_mat,
                                          annot_mat = annot_mat)
  
  ### Separate input matrices into each cell type
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegSepMat.R"))
  res.sep.mat <- cellCellRegSepMat(exp_mat = exp_mat,
                                   fet_mat = res.make.features$fet_mat_unnorm,
                                   annot_score_mat = res.make.features$annot_score_mat,
                                   HVG_opt = HVG_opt,
                                   HVG_extract_num = HVG_extract_num)
  
  ### Estimate by PLS regression
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegEst.R"))
  res.estimate <- cellCellRegEst(fet_mat_orig = res.make.features$fet_mat_unnorm,
                                 coord_mat = coord_mat,
                                 annot_score_mat = res.make.features$annot_score_mat,
                                 exp_mat_sep_list = res.sep.mat$exp_mat_sep_list,
                                 fet_mat_sep_list = res.sep.mat$fet_mat_sep_list,
                                 estimate_cell_type_opt = "all",
                                 cell_type_list = res.sep.mat$cell_type_list,
                                 dev_opt = dev_opt)
  
  ### Select significant variable
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegSelVar.R"))
  res.sel.var <- cellCellRegSelVar(res.estimate, res.sep.mat, HVG_extract_num,
                                   dev_opt = "kmeans")
  
  ### Get heatmap object
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegHeatmap.R"))
  res.heatmap <- cellCellRegHeatmap(res.sel.var, output_dir)

  ### Get bipartite graph object
  source(paste0(my_working_dir, "/v0.1.8/R/cellCellRegGraph.R"))
  res.graph <- cellCellRegGraph(res.sel.var, output_dir)
  
  
  return(list(res.make.features = res.make.features,
              res.sep.mat = res.sep.mat,
              res.estimate = res.estimate,
              res.sel.var = res.sel.var,
              res.heatmap = res.heatmap,
              res.graph = res.graph))
  
}
