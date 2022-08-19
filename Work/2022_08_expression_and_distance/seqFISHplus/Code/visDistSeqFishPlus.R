visDistSeqFishPlus <- function(res.sep.mat, tg = "Mag",
                               dc = c("Olig", "OPC", "astrocytes"),
                               dc_type = "1",
                               coord_mat, annot_mat){
  
  library(dplyr)
  library(ggplot2)
  
  opc_exp_mat <- res.sep.mat$exp_mat_sep_list[[12]]
  opc_tg_exp <- opc_exp_mat[, tg]
  
  if (sum(dc %in% "OPC") == 0){
    dc_2 <- append(dc, "OPC")
  } else {
    dc_2 <- dc
  }
  dc_2_flag <- annot_mat$cell_type %in% dc_2
  
  annot_mat_2 <- annot_mat[dc_2_flag,]
  coord_mat_2 <- coord_mat[dc_2_flag,]
  
  dist_mat <- as.matrix(stats::dist(coord_mat_2[, c("X", "Y")]))
  
  opc_flag <- annot_mat$cell_type %in% "OPC"
  opc_row <- rownames(annot_mat[opc_flag, ])
  
  dc_flag <- annot_mat$cell_type %in% dc
  dc_col <- rownames(annot_mat[dc_flag, ])
  
  dist_mat_2 <- dist_mat[opc_row, dc_col]
  dist_mat_2[dist_mat_2 == 0] <- 10000
  min_dist <- apply(dist_mat_2, 1, min)
  
  # Transformation to μm scale
  min_dist_2 <- min_dist * 0.106
  
  df <- as_tibble(cbind(opc_tg_exp, min_dist_2))
  colnames(df) <- c("Gene_expression", "Distance")
  
  if (dc_type == 1){
    xlab_txt <- c("Nearest distance to any of astrocytes, Olig and OPCs")
  } else if (dc_type == 2){
    xlab_txt <- c("Nearest distance to any of astrocytes and Olig")
  } else if (dc_type == 3) {
    xlab_txt <- "Nearest distance to OPCs"
  } else if (dc_type == 4) {
    xlab_txt <- "Nearest distance to astrocytes"
  } else if (dc_type == 5) {
    xlab_txt <- "Nearest distance to Olig"
  }
  
  p <- ggplot(df, aes(x = Distance, y = Gene_expression)) +
              geom_point() +
              theme_bw() +
              # ggtitle(tg) +
              xlim(0, 400) +
              ylab(paste0(tg, " expression in OPCs")) +
              xlab(xlab_txt)

  ggsave(paste0("~/CCPLS_paper/Work/2022_08_expression_and_distance/seqFISHplus/PDF/",
                tg, "_", dc_type, ".pdf"),
         p, width = 6, height = 3)
  ggsave(paste0("~/CCPLS_paper/Work/2022_08_expression_and_distance/seqFISHplus/PNG/",
                tg, "_", dc_type,  ".png"),
         p, width = 6, height = 3)
  
  
}

work_dir <- "~/CCPLS_paper/Work/2022_08_expression_and_distance/seqFISHplus/Code/"
coord_mat <- readRDS(paste0(work_dir, "coord_mat.rds"))
annot_mat <- readRDS(paste0(work_dir, "annot_mat.rds"))
res.sep.mat <- readRDS(paste0(work_dir, "res.sep.mat.rds"))

# Tmem98 / Abl1 / Dock6 / Mag / Kndc1
tg_list <- c("Tmem98", "Abl1", "Dock6", "Mag", "Kndc1")

# any of astrocytes, Olig, OPC
visDistSeqFishPlus(res.sep.mat, tg = "Tmem98",
                   dc = c("Olig", "OPC", "astrocytes"),
                   dc_type = "1",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Abl1",
                   dc = c("Olig", "OPC", "astrocytes"),
                   dc_type = "1",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Dock6",
                   dc = c("Olig", "OPC", "astrocytes"),
                   dc_type = "1",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Mag",
                   dc = c("Olig", "OPC", "astrocytes"),
                   dc_type = "1",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Kndc1",
                   dc = c("Olig", "OPC", "astrocytes"),
                   dc_type = "1",
                   coord_mat, annot_mat)

# any of astrocytes, Olig
visDistSeqFishPlus(res.sep.mat, tg = "Tmem98",
                   dc = c("Olig", "astrocytes"),
                   dc_type = "2",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Abl1",
                   dc = c("Olig", "astrocytes"),
                   dc_type = "2",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Dock6",
                   dc = c("Olig", "astrocytes"),
                   dc_type = "2",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Mag",
                   dc = c("Olig", "astrocytes"),
                   dc_type = "2",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Kndc1",
                   dc = c("Olig", "astrocytes"),
                   dc_type = "2",
                   coord_mat, annot_mat)

# OPC
visDistSeqFishPlus(res.sep.mat, tg = "Tmem98",
                   dc = c("OPC"),
                   dc_type = "3",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Abl1",
                   dc = c("OPC"),
                   dc_type = "3",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Dock6",
                   dc = c("OPC"),
                   dc_type = "3",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Mag",
                   dc = c("OPC"),
                   dc_type = "3",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Kndc1",
                   dc = c("OPC"),
                   dc_type = "3",
                   coord_mat, annot_mat)

# astrocytes
visDistSeqFishPlus(res.sep.mat, tg = "Tmem98",
                   dc = c("astrocytes"),
                   dc_type = "4",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Abl1",
                   dc = c("astrocytes"),
                   dc_type = "4",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Dock6",
                   dc = c("astrocytes"),
                   dc_type = "4",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Mag",
                   dc = c("astrocytes"),
                   dc_type = "4",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Kndc1",
                   dc = c("astrocytes"),
                   dc_type = "4",
                   coord_mat, annot_mat)

# Olig
visDistSeqFishPlus(res.sep.mat, tg = "Tmem98",
                   dc = c("Olig"),
                   dc_type = "5",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Abl1",
                   dc = c("Olig"),
                   dc_type = "5",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Dock6",
                   dc = c("Olig"),
                   dc_type = "5",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Mag",
                   dc = c("Olig"),
                   dc_type = "5",
                   coord_mat, annot_mat)
visDistSeqFishPlus(res.sep.mat, tg = "Kndc1",
                   dc = c("Olig"),
                   dc_type = "5",
                   coord_mat, annot_mat)