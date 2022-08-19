visSpatSeqScope <- function(res.sep.mat, tg = "Gpx1",
                               coord_mat, annot_mat){
  
  library(dplyr)
  library(ggplot2)
  
  imm_exp_mat <- res.sep.mat$exp_mat_sep_list[[3]]
  imm_tg_exp <- imm_exp_mat[, tg]
  plot_list <- c("B.cell-IgA", "B.cell-Immature")
  plot_flag <- annot_mat[, 3] %in% plot_list
  annot_mat_2 <- annot_mat
  annot_mat_2[!plot_flag, 3] <- "Others"
  
  imm_flag <- annot_mat_2[, 3] == "B.cell-Immature"
  imm_coord <- coord_mat[imm_flag, ]
  imm_exp_coord <- cbind(as.double(imm_tg_exp), imm_coord) %>% as_tibble()
  
  iga_flag <- annot_mat_2[, 3] %in% "B.cell-IgA"
  iga_coord <- as_tibble(coord_mat[iga_flag, ])
  iga_coord_2 <- cbind(iga_coord, rep("B cell-IgA", nrow(iga_coord)))
  colnames(iga_coord_2) <- c("cell_ID", "X", "Y", "Cell_type")
  
  oth_coord <- as_tibble(coord_mat[!plot_flag, ])
  oth_coord_2 <- cbind(oth_coord, rep("Others", nrow(oth_coord)))
  colnames(oth_coord_2) <- c("cell_ID", "X", "Y", "Cell_type")
  
  p <- ggplot() +
    geom_point(data = imm_exp_coord, aes(x = x_coord, y = y_coord, color = imm_tg_exp)) +
    scale_color_viridis_c() +
    geom_point(data = iga_coord_2, aes(x = X, y = Y, fill = Cell_type),
               color = "black", shape = "triangle") +
    geom_point(data = oth_coord_2, aes(x = X, y = Y),
               color = "grey", alpha = 0.1, shape = "plus") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    ggtitle(paste0("Spatial distribution of ", tg, " expression in B cell-Immature")) +
    xlab("X coordinate") +
    ylab("Y coordinate")
  
  ggsave(paste0("~/CCPLS_paper/Work/2022_08_spatial_visualization/Seq-Scope/PDF/", tg, ".pdf"),
         p, width = 14, height = 7)
  ggsave(paste0("~/CCPLS_paper/Work/2022_08_spatial_visualization/Seq-Scope/PNG/", tg, ".png"),
         p, width = 14, height = 7)
  
}

work_dir <- "~/CCPLS_paper/Work/2022_08_spatial_visualization/Seq-Scope/Code/"
coord_mat <- readRDS(paste0(work_dir, "coord_mat.rds"))
annot_mat <- readRDS(paste0(work_dir, "annot_mat.rds"))
res.sep.mat <- readRDS(paste0(work_dir, "res.sep.mat.rds"))
res.sel.var <- readRDS(paste0(work_dir, "res.sel.var.rds"))

gc <- res.sel.var$gene_cluster_vec_list[[3]]
tg_list <- names(gc[gc == 4])

# Tbc1d20 / Hnf4a / Sox18 / Rapgef2 / Nup210l / Myadm / Sidt2 / Gpx1 /
# Psap / Pcdh15 / Xbp1 / Foxa1 / Atf4 / Fshr / Yipf6
# tg_list <- c("Tbc1d20", "Hnf4a", "Sox18", "Rapgef2", "Nup210l", "Myadm", "Sidt2", "Gpx1",
#              "Psap", "Pcdh15", "Xbp1", "Foxa1", "Atf4", "Fshr", "Yipf6")
for (i in 1:length(tg_list)){
  visSpatSeqScope(res.sep.mat, tg = tg_list[i],
                  coord_mat, annot_mat)
}