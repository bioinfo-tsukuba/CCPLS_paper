visSpatSeqFishPlus <- function(res.sep.mat, tg = "Mag",
                               coord_mat, annot_mat){
  
  library(dplyr)
  library(ggplot2)
  
  opc_exp_mat <- res.sep.mat$exp_mat_sep_list[[12]]
  opc_tg_exp <- opc_exp_mat[, tg]
  plot_list <- c("Olig", "OPC", "astrocytes")
  plot_flag <- annot_mat$cell_type %in% plot_list
  annot_mat_2 <- annot_mat
  annot_mat_2$cell_type[!plot_flag] <- "Others"
  
  opc_flag <- annot_mat_2$cell_type == "OPC"
  opc_coord <- coord_mat[opc_flag, ]
  opc_exp_coord <- as_tibble(cbind(opc_tg_exp, opc_coord))
  
  olig_flag <- annot_mat_2$cell_type %in% "Olig"
  olig_coord <- as_tibble(coord_mat[olig_flag, ])
  olig_coord_2 <- cbind(olig_coord, rep("Olig", nrow(olig_coord)))
  colnames(olig_coord_2) <- c("cell_ID", "X", "Y", "Cell_type")
  
  ast_flag <- annot_mat_2$cell_type %in% "astrocytes"
  ast_coord <- as_tibble(coord_mat[ast_flag, ])
  ast_coord_2 <- cbind(ast_coord, rep("astrocytes", nrow(ast_coord)))
  colnames(ast_coord_2) <- c("cell_ID", "X", "Y", "Cell_type")
  
  oth_coord <- as_tibble(coord_mat[!plot_flag, ])
  oth_coord_2 <- cbind(oth_coord, rep("Others", nrow(oth_coord)))
  colnames(oth_coord_2) <- c("cell_ID", "X", "Y", "Cell_type")
  
  p <- ggplot() +
    geom_point(data = opc_exp_coord, aes(x = X, y = Y, color = opc_tg_exp)) +
    scale_color_viridis_c() +
    geom_point(data = olig_coord_2, aes(x = X, y = Y, fill = Cell_type),
               alpha = 0.6, color = "black", shape = "triangle") +
    geom_point(data = ast_coord_2, aes(x = X, y = Y, fill = Cell_type),
               alpha = 0.6, color = "black", shape = "square") +
    geom_point(data = oth_coord_2, aes(x = X, y = Y, fill = Cell_type),
               alpha = 0.1, shape = "plus") +
    theme_bw() +
    theme(aspect.ratio = 0.5) +
    ggtitle(paste0("Spatial distribution of ", tg, " expression in OPC")) +
    xlab("X coordinate") +
    ylab("Y coordinate")
  
  ggsave(paste0("~/CCPLS_paper/Work/2022_08_spatial_visualization/seqFISHplus/", tg, ".pdf"),
         p, width = 14, height = 7)
  
  
}

work_dir <- "~/CCPLS_paper/Work/2022_08_spatial_visualization/seqFISHplus/Code/"
coord_mat <- readRDS(paste0(work_dir, "coord_mat.rds"))
annot_mat <- readRDS(paste0(work_dir, "annot_mat.rds"))
res.sep.mat <- readRDS(paste0(work_dir, "res.sep.mat.rds"))

# Tmem98 / Abl1 / Dock6 / Mag / Kndc1
tg_list <- c("Tmem98", "Abl1", "Dock6", "Mag", "Kndc1")
for (i in 1:length(tg_list)){
  visSpatSeqFishPlus(res.sep.mat, tg = tg_list[i],
                     coord_mat, annot_mat)
}