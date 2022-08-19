convertAnnotMatV13.1 <- function(annot_mat_seqFISHplus){

  annot_mat_pre <- annot_mat_seqFISHplus

  # A
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "L5 eNeuron", "cell_type"] <- "A"
  
  # B
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "L6 eNeuron", "cell_type"] <- "B"

  # C
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "Olig", "cell_type"] <- "C"
  
  # D
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "Adarb2 iNeuron", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "mural", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "L4 eNeuron", "cell_type"] <- "D"

  annot_mat_pre[annot_mat_pre[,"cell_type"] == "L2/3 eNeuron", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "Lhx6 iNeuron", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "microglia", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "endothelial", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "astrocytes", "cell_type"] <- "D"
  annot_mat_pre[annot_mat_pre[,"cell_type"] == "OPC", "cell_type"] <- "D"
  
  annot_mat <- annot_mat_pre
  
  return(annot_mat)
  
}