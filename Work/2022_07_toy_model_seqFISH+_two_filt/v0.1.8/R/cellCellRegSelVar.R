sigCoef <- function(res.pls,
                    exp_mat_norm,
                    fet_mat_norm,
                    p_value_thresh = 0.05){
  
  ## Prepare variable
  res_pls <- res.pls$res_pls
  opt_comp_num <- res.pls$opt_comp_num

  x_score_mat <- res_pls$scores
  y_score_mat <- res_pls$Yscores
  
  ## Calculate sig_coef_mat
  
  # Declare sig_coef_mat
  sig_coef_mat <- matrix(0, nrow = nrow(res_pls$coefficient[,,1]),
                         ncol = ncol(res_pls$coefficients[,,1]))
  rownames(sig_coef_mat) <- rownames(res_pls$coefficient[,,1])
  colnames(sig_coef_mat) <- colnames(res_pls$coefficient[,,1])
  sig_coef_mat_raw <- sig_coef_mat
  
  # Add each component
  for (calc_comp_ind in 1:opt_comp_num){
    
    if (calc_comp_ind > 1){
      sig_coef_mat_comp_raw <- res_pls$coefficients[,,calc_comp_ind] - res_pls$coefficients[,,calc_comp_ind - 1]
    } else {
      sig_coef_mat_comp_raw <- res_pls$coefficients[,,calc_comp_ind]
    }
    
    # Extract non-significant features
    del_fet_ind_vec <- returnDelFetInd(fet_mat_norm = fet_mat_norm,
                                       x_score_mat = x_score_mat,
                                       calc_comp_num = calc_comp_ind,
                                       p_value_thresh = p_value_thresh)
    
    # Extract non-significant genes
    del_gene_ind_vec <- returnDelGeneInd(exp_mat_norm = exp_mat_norm,
                                         y_score_mat = y_score_mat,
                                         calc_comp_num = calc_comp_ind,
                                         p_value_thresh = p_value_thresh)
    
    # Replace non-significant value with 0
    sig_coef_mat_comp <- sig_coef_mat_comp_raw
    sig_coef_mat_comp[del_fet_ind_vec, ] <- 0
    sig_coef_mat_comp[, del_gene_ind_vec] <- 0
    
    sig_coef_mat_raw <- sig_coef_mat_raw + sig_coef_mat_comp_raw
    sig_coef_mat <- sig_coef_mat + sig_coef_mat_comp
    
  }
  
  return(list(sig_coef_mat = sig_coef_mat,
              sig_coef_mat_raw = sig_coef_mat_raw))
  
}

returnDelFetInd <- function(fet_mat_norm,
                            x_score_mat,
                            calc_comp_num,
                            p_value_thresh){
  
  # Remove self
  self_col_ind <- grep("self", colnames(fet_mat_norm))
  fet_mat_norm2 <- fet_mat_norm[, -self_col_ind]
  
  fet_p_vec <- c()
  del_fet_ind_vec <- c()
  
  fet_num <- ncol(fet_mat_norm2)
  
  for (fet_ind in 1:fet_num){
    fet_p_value <- cor.test(fet_mat_norm2[, fet_ind], x_score_mat[, calc_comp_num])$p.value
    fet_p_vec <- append(fet_p_vec, fet_p_value)
  }
  
  fet_q_vec <- p.adjust(fet_p_vec, method = "BH")
  del_fet_bin <- rep(1, fet_num) * (fet_q_vec >= p_value_thresh)
  
  for (fet_ind in 1:fet_num){
    if (del_fet_bin[fet_ind] == 1)
    del_fet_ind_vec <- c(del_fet_ind_vec, fet_ind)
  }
  
  return(del_fet_ind_vec) 
  
}

returnDelGeneInd <- function(exp_mat_norm,
                             y_score_mat,
                             calc_comp_num,
                             p_value_thresh){
  gene_p_vec <- c()
  del_gene_ind_vec <- c()
  
  gene_num <- ncol(exp_mat_norm)
  
  for (gene_ind in 1:gene_num){
    gene_p_value <- cor.test(exp_mat_norm[, gene_ind], y_score_mat[, calc_comp_num])$p.value
    gene_p_vec <- append(gene_p_vec, gene_p_value)
  }
  
  gene_q_vec <- p.adjust(gene_p_vec, method = "BH")
  del_gene_bin <- rep(1, gene_num) * (gene_q_vec >= p_value_thresh)
  
  for (gene_ind in 1:gene_num){
    if (del_gene_bin[gene_ind] == 1)
      del_gene_ind_vec <- c(del_gene_ind_vec, gene_ind)
  }

  return(del_gene_ind_vec) 
  
}

cellCellRegSelVar <- function(res.estimate, res.sep.mat, HVG_extract_num, dev_opt){
  
  set.seed(123)
  
  cell_type_num <- length(res.estimate$cell_type_list)
  
  res_xmeans_list <- vector("list", length = cell_type_num)
  
  sig_coef_mat_non_zero_list <- vector("list", length = cell_type_num)
  sig_coef_mat_with_zero_list <- vector("list", length = cell_type_num)
  
  sig_coef_mat_list <- vector("list", length = cell_type_num)
  ave_cluster_coef_mat_list <- vector("list", length = cell_type_num)
  
  sig_coef_mat_bin_2_list <- vector("list", length = cell_type_num)
  gene_cluster_vec_list <- vector("list", length = cell_type_num)
  
  all_zero_flag_list <- vector("list", length = cell_type_num)
  
  for (cell_type_ind in 1:cell_type_num){
    
    # Judge model was built or not.
    if (res.estimate$CCPLS_result[[cell_type_ind]][[1]] != "No model was built."){
      
      # Get significant coefficient
      res.sig.coef <- sigCoef(res.pls = res.estimate$CCPLS_result[[cell_type_ind]][[2]],
                             exp_mat_norm = scale(res.sep.mat$exp_mat_sep_list[[cell_type_ind]]),
                             fet_mat_norm = res.sep.mat$fet_mat_sep_list[[cell_type_ind]])
      
      sig_coef_mat <- res.sig.coef$sig_coef_mat
      sig_coef_mat_raw <- res.sig.coef$sig_coef_mat_raw
      
      non_sig_gene_col <- colnames(sig_coef_mat)[apply(sig_coef_mat, 2, sum) == 0]
      sig_gene_col <- colnames(sig_coef_mat)[apply(sig_coef_mat, 2, sum) != 0]
      non_sig_fet_row <- rownames(sig_coef_mat)[apply(sig_coef_mat, 1, sum) == 0]
      sig_fet_row <- rownames(sig_coef_mat)[apply(sig_coef_mat, 1, sum) != 0]

      sig_coef_mat_bin <- matrix(0, nrow(sig_coef_mat), ncol(sig_coef_mat))
      rownames(sig_coef_mat_bin) <- rownames(sig_coef_mat)
      colnames(sig_coef_mat_bin) <- colnames(sig_coef_mat)
      

      saveRDS(sig_coef_mat_raw, paste0(my_working_dir, output_dir, "/sig_coef_mat_raw.rds"))
      saveRDS(sig_coef_mat, paste0(my_working_dir, output_dir, "/sig_coef_mat_step_1.rds"))
      
      for (fet_id in 1:nrow(sig_coef_mat)){
          
        if (length(non_sig_gene_col) != 0){
          null_dist_fet <- ecdf(sig_coef_mat_raw[fet_id, non_sig_gene_col])
          fet_p_vec <- 1 - null_dist_fet(sig_coef_mat_raw[fet_id,])
          names(fet_p_vec) <- colnames(sig_coef_mat_raw)
            
          fet_q_vec <- p.adjust(fet_p_vec, method = "BH") 
            
          sig_coef_vec_fet <- sig_coef_mat[fet_id,]
          sig_plus_flag <- fet_q_vec < 0.05 & sig_coef_vec_fet > 0
          sig_minus_flag <- fet_q_vec < 0.05 & sig_coef_vec_fet < 0
          sig_coef_mat_bin[fet_id, sig_plus_flag] <- 1
          sig_coef_mat_bin[fet_id, sig_minus_flag] <- -1
          
          # Check removed coefficients
          # sig_coef_mat_raw - sig_coef_matがfiltering 1
          fo_flag <- fet_q_vec >= 0.05 & sig_coef_vec_fet != 0
          fo_coef_vec <- sig_coef_mat[fet_id, fo_flag]
          saveRDS(fo_coef_vec, paste0(my_working_dir, output_dir, "/fo_coef_vec_fet_",
                                      rownames(sig_coef_mat)[fet_id],".rds"))
          
        }
          
      }
        
      sig_coef_mat_bin_2 <- sig_coef_mat_bin
      sig_coef_mat_bin_2[, non_sig_gene_col] <- 0
      sig_coef_mat_bin_2[non_sig_fet_row, ] <- 0
      sig_coef_mat_bin_2_non_zero <- sig_coef_mat_bin_2[, colSums(sig_coef_mat_bin_2) != 0]
      
      sig_coef_mat[sig_coef_mat_bin_2 == 0] <- 0
      saveRDS(sig_coef_mat, paste0(my_working_dir, output_dir, "/sig_coef_mat_step_2.rds"))
      
      sig_coef_mat_with_zero_list[[cell_type_ind]] <- sig_coef_mat
      sig_coef_mat_non_zero <- sig_coef_mat[, colSums(sig_coef_mat_bin_2) != 0]
      # zero_gene_num <- ncol(sig_coef_mat) - ncol(sig_coef_mat_non_zero)
      
      sig_coef_mat_non_zero_list[[cell_type_ind]] <- sig_coef_mat_non_zero
      
      all_zero_flag <- sum(sig_coef_mat_bin_2) == 0
      all_zero_flag_list[[cell_type_ind]] <- all_zero_flag
      
      if (!all_zero_flag){
        
        if (dev_opt == "xmeans"){
          
          res.xmeans <- xmeans(t(x = sig_coef_mat_non_zero_list[[cell_type_ind]]))
          cluster_num <- max(res.xmeans$cluster)
          gene_cluster_vec <- res.xmeans$cluster
          res_xmeans_list[[cell_type_ind]] <- res.xmeans
          gene_cluster_vec_list[[cell_type_ind]] <- gene_cluster_vec
        
        } else if (dev_opt == "kmeans"){
          
          set.seed(1)
          
          mat <- scale(t(sig_coef_mat_non_zero[rowSums(sig_coef_mat_non_zero) != 0,]))
          
          k.max <- 15
          
          if ((nrow(mat) - 1) < 15){
            k.max <- nrow(mat) - 1
          }
          
          sil_width <- purrr::map_dbl(2:k.max,  function(k){
            model <- cluster::pam(x = mat, k = k)
            model$silinfo$avg.width
          })
          
          cluster_num <- which.max(sil_width) + 1
          res.kmeans <- kmeans(mat, cluster_num)
          gene_cluster_vec <- res.kmeans$cluster
          gene_cluster_vec_list[[cell_type_ind]] <- gene_cluster_vec
          
        } else if (dev_opt == "binary"){

          gene_bin_vec_list <- vector("list", ncol(sig_coef_mat_bin_2_non_zero))
          for (gene_id in 1:ncol(sig_coef_mat_bin_2_non_zero)){
            gene_bin_vec_list[[gene_id]] <- as.character(sig_coef_mat_bin_2_non_zero[,gene_id])
          }
          
          cluster_list <- unique(gene_bin_vec_list)
          cluster_num <- length(cluster_list)
          
          gene_cluster_vec <- c()
          for (cluster_id in 1:cluster_num){
            for (gene_id in 1:length(gene_bin_vec_list)){
              judge_flag <- nrow(sig_coef_mat_bin_2) == sum(gene_bin_vec_list[[gene_id]] == cluster_list[[cluster_id]])
              if (judge_flag){
                gene_cluster_vec <- append(gene_cluster_vec, cluster_id)
              }
            }
          }
          
          names(gene_cluster_vec) <- colnames(sig_coef_mat_bin_2_non_zero)
          gene_cluster_vec_list[[cell_type_ind]] <- gene_cluster_vec
          
        }

        cluster_name <- c()
        ave_cluster_coef_mat_non_zero <- matrix(0, nrow = nrow(sig_coef_mat_non_zero_list[[cell_type_ind]]),
                                                ncol = cluster_num)
        
        for (cluster_ind in 1:cluster_num){
          
          cluster_gene <- names(gene_cluster_vec[gene_cluster_vec == cluster_ind])
          cluster_coef_mat_non_zero <- sig_coef_mat_non_zero_list[[cell_type_ind]][, cluster_gene]
          
          if (!is.null(ncol(cluster_coef_mat_non_zero))){ # judge column is 1 or not
            ave_cluster_coef_mat_non_zero[, cluster_ind] <- apply(cluster_coef_mat_non_zero,
                                                                  1, mean)
          } else {
            ave_cluster_coef_mat_non_zero[, cluster_ind] <- cluster_coef_mat_non_zero
          }
          
          cluster_name <- append(cluster_name, paste0("SRGs Cluster ", cluster_ind))
          
        }
        
        rownames(ave_cluster_coef_mat_non_zero) <- rownames(sig_coef_mat_non_zero)
        colnames(ave_cluster_coef_mat_non_zero) <- cluster_name
        
        if (ncol(sig_coef_mat_non_zero) <= HVG_extract_num){
          non_DEG_vec <- matrix(0, nrow = nrow(ave_cluster_coef_mat_non_zero), ncol = 1)
          colnames(non_DEG_vec) <- "non-SRGs"
          ave_cluster_coef_mat <- cbind(ave_cluster_coef_mat_non_zero,
                                        non_DEG_vec)
          ave_cluster_coef_mat_list[[cell_type_ind]] <- ave_cluster_coef_mat
        }

        sig_coef_mat_non_zero_ordered <- sig_coef_mat_non_zero_list[[cell_type_ind]][, order(gene_cluster_vec)]
        
        # Remove zero feature
        sig_coef_mat <- sig_coef_mat_non_zero_ordered[rowSums(sig_coef_mat_non_zero_ordered) != 0, ]
        
        row_name_rm <- c()
        for (row_ind in 1:nrow(sig_coef_mat)){
          row_name_rm <- append(row_name_rm, strsplit(rownames(sig_coef_mat)[row_ind], split = "neig_")[[1]][2])
        }
        rownames(sig_coef_mat) <- row_name_rm
        
        sig_coef_mat_list[[cell_type_ind]] <- sig_coef_mat
        sig_coef_mat_bin_2_list[[cell_type_ind]] <- sig_coef_mat_bin_2
      
      } else { # in case of all zero
        
        res_xmeans_list [[cell_type_ind]] <- "NULL"
        gene_cluster_vec_list[[cell_type_ind]] <- "NULL"
        ave_cluster_coef_mat_list[[cell_type_ind]] <- "NULL"
        sig_coef_mat_list[[cell_type_ind]] <- "NULL"
        sig_coef_mat_bin_2_list[[cell_type_ind]] <- "NULL"

      }
      
    }

  }
    
  return(list(res_xmeans_list = res_xmeans_list,
              ave_cluster_coef_mat_list = ave_cluster_coef_mat_list,
              sig_coef_mat_list = sig_coef_mat_list,
              sig_coef_mat_with_zero_list = sig_coef_mat_with_zero_list,
              sig_coef_mat_bin_2_list = sig_coef_mat_bin_2_list,
              gene_cluster_vec_list = gene_cluster_vec_list,
              all_zero_flag_list = all_zero_flag_list,
              cell_type_list = res.estimate$cell_type_list,
              cell_type_model_SG = res.estimate$cell_type_model_SG))
  
}

# $Id: xmeans.prog,v 1.23 2012/04/12 11:21:12 tunenori Exp tunenori $
#
# X-MEANS Clustering
#
# Description:
#
#      Perform x-maens non-hierarchal clustering on a data matrix.
# 
# Usage:
# 
#      xmeans(x, ik = 2, iter.max = 10, pr.proc = F, 
#		ignore.covar = T, merge.cls = F)
# 
# Arguments:
#	x: A numeric matrix of data, or an object that can be coerced to
#          such a matrix (such as a numeric vector or a data frame with
#          all numeric columns).
#
#	ik: The initial  number of clusters applied to kmeans().
#	   As xmeans calls kmeans recursively, `ik' should be sufficient
#	   small.
#
# 	iter.max: The maximum of iterations allowed.
#
#	pr.proc: logical: If 'TRUE' the system outputs the processing status.
#
# 	ignore.covar: logical: If 'TRUE', covariances of cluster data are
#	   ignored. For saving of the time, 'TRUE' is set as the defalut.
#
#	merge.cls: logical: If 'TRUE', some clusters may be merged into another
#	   clusters after iterative division. 
#
# Value:
#
#     An object of class 'xmeans' which is a list with components:
#
#	cluster: A vector of integers indicating the cluster to which each
#          point is allocated.
#
#	centers: A matrix of cluster centres.
#
#	size: The number of points in each cluster. When `merge.cls' is TRUE,
#	   some elements may be zero.
#
# References:
#
#	Ishioka, T. (2005): ``An Expansion of X-means for Automatically
#	Determining the Optimal Number of Clusters," The Fourth IASTED
#	International Conference on Computational Intelligence (CI 2005),
#	Calgary Canada, July 4-6, pp.91-96. 
#	http://www.rd.dnc.ac.jp/%7Etunenori/doc/487-053.pdf
#
#	Ishioka, T. (2000): ``Extended K-means with an Efficient Estimation
#	of the number of Clusters,'' Intelligent Data Engineering and
#	Automated Learning --- IDEAL 2000, Second International Conference,
#	Shatin, N.T., Hong Kong, China, December 2000, proceedings 17--22.
#	(Lecture Notes in Computer Science 1983, Kwong Sak Leung, Lai-Wan
#	Chan, Helen Meng (Eds.), Springer, 17--22, 2000) 
#	http://www.rd.dnc.ac.jp/%7Etunenori/doc/xmeans_ideal2000.pdf
#
# Examples:
#
#	xmeans(iris[,-5],  merge.cls=T)
#	plot(cmdscale(dist(iris[,-5])), cex=2, pch=as.numeric(iris[,5]))
#

xmeans <- function(x, ik = 2, iter.max = 10, pr.proc = F, ignore.covar = T, merge.cls = F){
  if (ik < 2) 
    ik <- 2
  x <- as.matrix(x)
  p <- ncol(x) # p-dimensional multivariate
  if (ignore.covar){
    q <- 2 * p # number of parameters; mean and var for each "p"
  }else{
    q <- p * (p+3) / 2	# integer
  }
  cl<- kmeans(x,ik,iter.max)
  cl.sub <- list()
  
  for (i in 1:ik){ # for each cluster
    y.ok <- (cl$cluster == i) 	# i-th cluster or not
    yi <- matrix(x[y.ok], ncol=p) 	# extract i-th cluster
    zi <- yi 	# save the data for graphics
    yi.centers <- cl$centers[i,]
    zi.centers <- yi.centers
    yi.cluster <- cl$cluster[(cl$cluster == i)]
    yi.cluster <- rep(1, length(yi.cluster)) 
    # sub-cluster number should begin from 1
    
    k1 <- 1		# cluster number
    k2 <- k1 + 1
    bic.prior <- NULL
    stack <- list()	# divided and unproceeded data are stacked
    lnL0 <- lnL(yi, yi.centers, ignore.covar)
    yi.lnL <- lnL0$lnL
    yi.detVx <- lnL0$detVx
    
    repeat{
      
      # go through at least 1 time; 
      # y$subcluster exist...
      if (pr.proc)	cat (paste("k1 =", k1, ", k2 =", k2,"\n"))
      if (nrow(yi) == 1){ # sample size is 1
        break
      }
      y <- split2cls(yi, yi.centers, q, bic.prior, lnL.prior, detVx.prior, iter.max, ignore.covar)
      if (y$continue){ # splitting continue 
        yi.cluster <-
          updtCrusterNum(y$continue, yi.cluster, k1, k2, y$subcluster)
        zi.centers <-
          updtCenters(y$continue, zi.centers, k1, k2, y$centers)
        yi.lnL <-
          updtlnL(y$continue, yi.lnL, k1, k2, y$lnL.post)
        yi.detVx <-
          updtdetVx(y$continue, yi.detVx, k1, k2, y$detVx.post)
      }
      
      if (pr.proc) print(y$subcluster)
      if (pr.proc){ print(y$bic.prior)
        print(y$bic.post)
        # print(y$lnL.prior)	# for debug
        # print(y$lnL.post)	# for debug
        print(y$continue) }
      # cat("zi.centers=\n")	# for debug
      # print(zi.centers)	# for debug
      if (!y$continue){	# no-node
        if ((nstack <- length(stack))){ # there are stacked data
          # extract the stacked data
          if (pr.proc)
            cat(paste("extract the stacked data (", nstack, ")...\n"))
          yi <- stack[[nstack]]$data
          yi.centers <- stack[[nstack]]$centers
          bic.prior <- stack[[nstack]]$bic
          lnL.prior <- stack[[nstack]]$lnL
          detVx.prior <- stack[[nstack]]$detVx
          k1 <- stack[[nstack]]$cls
          k2 <- k2 # unchanged
          # delete the data set
          if (nstack > 1){
            stack <- stack[1:(nstack-1)]
          }else{
            stack <- list() # no stacked data
          }
          next;
        }
        # no node and no stack
        if (pr.proc)	cat ("no node and no stack...\n")
        break;
      }
      # splitting continues...
      y1 <- y$clj1	# data
      y2 <- y$clj2
      yi.ctr1 <- y$centers[1,]	# centers
      yi.ctr2 <- y$centers[2,]
      bic.prior1 <- y$bic.post[1]	# bic
      bic.prior2 <- y$bic.post[2]
      lnL.prior1 <- y$lnL.post[1]	# lnL
      lnL.prior2 <- y$lnL.post[2]
      detVx.prior1 <- y$detVx.post[1]	# detVx
      detVx.prior2 <- y$detVx.post[2]
      
      # one-hand repeats recursively...
      yi <- y1
      yi.centers <- yi.ctr1
      bic.prior <- bic.prior1
      lnL.prior <- lnL.prior1
      detVx.prior <- detVx.prior1
      # other-hand is stacked... 
      if (pr.proc)	cat ("stacking ...\n")
      stack <- c(stack,
                 list(list(data=y2, centers=yi.ctr2,
                           bic=bic.prior2, lnL=lnL.prior2, detVx=detVx.prior2, cls=k2)))
      # inclement the cluster number 
      k2 <- k2 + 1
      
    } # end of repeat
    
    # splitting done ...
    if (pr.proc){
      cat ("splitting done...\n")
      cat (paste("main cluster =",i,"*******\n"))
    }
    cl.sub <- c(cl.sub, list(list(cluster = yi.cluster,
                                  centers = zi.centers, lnL = yi.lnL, detVx = yi.detVx,
                                  size = tabulate(yi.cluster))))
    if (pr.proc){
      print(cl.sub[[i]])
      plot(zi, col=yi.cluster)
      if (is.vector(zi.centers))
        points(zi.centers[1], zi.centers[2], pch=8)
      else # array
        points(zi.centers,col=1:(length(zi.centers)/p),pch=8)
    }
  }
  if (pr.proc)	print(cl.sub)
  xcl <- mergeResult(cl, cl.sub, ik)
  
  if (merge.cls == F) {
    return(list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size))
  }
  
  # merge after progressive dividing
  #
  if (pr.proc) cat("merging after progressive dividing ...\n")
  
  k <- length(xcl$size)	# final cluster number
  if (k <= 2){	# minimum cluster number should be 2
    if (pr.proc) cat("merging skipped ...\n")
    return(list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size))
  }
  if (pr.proc){
    cat("xcl$detVx=")
    print(xcl$detVx)
    cat("xcl$size=")
    print(xcl$size)
  }
  
  klist <- sort.list(xcl$size) # "small" to "large" order of xcl$detVx list
  if (pr.proc) print(klist)
  for (i in 1:(k-1)){ 
    for (j in (i+1):k){ 
      k1 = klist[i]
      k2 = klist[j]
      if (pr.proc) cat(paste("inspecting the clusters", k1,"and", k2,"\n"))
      
      z <- mergedBIC(x, xcl, k1, k2, q, ignore.covar, pr.proc)
      if (z$ret == F){
        # k1 or k2 has been merged.
        # skip this roop
        if (pr.proc) cat("skipping... k1=", k1, "k2=", k2,"\n")
        next
      }
      if (z$bicdiv > z$bicmgd){
        # we prefer merged model.
        # replace larger cls. number to smaller cls. number
        if (pr.proc) cat("replace cls.", k2, "to", k1,"\n")
        xcl$cluster <- replace(xcl$cluster, (xcl$cluster == k2), k1)
        xcl$size[k1] <- xcl$size[k1] + xcl$size[k2]
        xcl$size[k2] <- 0
        xcl$lnL[k1] <- z$lnLmgd
        xcl$lnL[k2] <- 0
        xcl$detVx[k1] <- z$detVxmgd
        xcl$detVx[k2] <- 0
        xcl$centers[k1,] <- z$ctrmgd
        xcl$centers[k2,] <- 0
        
      }
    }
  }
  list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size)
}




# marge the result of sub-clustering;
# cluster numbers by first kmeans should be renumbered;
# the other centers and sizes are simply added.
# cl: the result of first kmeans
# cl.sub: the result of subclustering
# ik: cluster number adopted to kmeans.
mergeResult <- function(cl, cl.sub, ik){
  cluster <- cl$cluster	# main cluster
  centers <- NULL
  size <- NULL
  lnL <- NULL
  detVx <- NULL
  
  k <- 0	# uniq cluster numbers; k should be decremental. 
  for (i in 1:ik)
    k <- k + length(cl.sub[[i]]$size)
  kk <- k
  
  for (i in ik:1){	# loop for main clusters obtained by kmeans
    xsub <- cl.sub[[i]]$cluster
    iki <- ik -i +1 
    centers <- rbind(centers, cl.sub[[iki]]$centers)
    size <- c(size, cl.sub[[iki]]$size)
    lnL <- c(lnL, cl.sub[[iki]]$lnL)
    detVx <- c(detVx, cl.sub[[iki]]$detVx)
    
    for (j in length(cl.sub[[i]]$size):1){ # loop for subclusters
      xsub <- replace(xsub, (xsub == j), k)
      k <- k -1
    }
    cluster <- replace(cluster, (cluster == i), xsub)
  }
  if (k != 0) stop("mergeResult: assertion failed (k = 0)...")
  dimnames(centers) <- list(1:kk, NULL)
  list(cluster = cluster, centers = centers, lnL = lnL, detVx = detVx, size = size)
}


# update the cluster number by using the result of "split2cls()"
# continue: no splitting
# v: cluster numbers vector for initial cluster.
# k1: cluster numbers should be updated; "k1" becomes "k1" and "k2"
# xsub: sub-cluster numbers vector of "v" whose value is "k1";
#	given "xsub" have 1 or 2.
updtCrusterNum <- function(continue, v, k1, k2, xsub){
  if (!is.vector(v)) 
    return(xsub)
  if (!continue)
    return(v)
  if (k1 == k2)
    stop("updtCrusterNum() : k1 and k2 should differ.")
  
  # below is same algorithm; explicit array operation is slow in R.
  # j <- 1
  # for (i in 1:length(v)){
  #	if (v[i] == k1){
  #		if (xsub[j] == 2)
  #			v[i] <- k2
  #		j <- j + 1
  #	}
  # }
  # end of algorithm
  xsub <- replace(xsub, (xsub == 2), k2) # changed
  xsub <- replace(xsub, (xsub == 1), k1) # unchanged
  v <- replace(v, (v == k1), xsub)
}


# update the cluster centers by using the result of "split2cls()"
# continue: no update
# org.centers: original centers matrix
# divided.centers: divided centers matrix; it has 2 rows.
updtCenters <- function(continue, org.centers, k1, k2, divided.centers){
  if (!is.matrix(org.centers)) 
    return(divided.centers)
  if (!continue)
    return(org.centers)
  if (k1 == k2)
    stop("updtCenters() : k1 and k2 should differ.")
  
  z <- NULL
  for (i in 1:max(k2, nrow(org.centers))){
    if (i == k1)
      z <- rbind(z, divided.centers[1,])
    else if (i == k2)
      z <- rbind(z, divided.centers[2,])
    else
      z <- rbind(z, org.centers[i,])
  }
  z
}

# update the lnL by using the result of "split2cls()"
# continue: no update
# org.lnL: original lnL vector
# divided.lnL: divided lnL vector having 2 elements.
updtlnL <- function(continue, org.lnL, k1, k2, divided.lnL){
  if (!is.vector(org.lnL))
    return(divided.lnL)
  if (!continue)
    return(org.lnL)
  if (k1 == k2)
    stop("updtlnL() : k1 and k2 should differ.")
  
  z <- NULL
  for (i in 1:max(k2, length(org.lnL))){
    if (i == k1)
      z <- c(z, divided.lnL[1])
    else if (i == k2)
      z <- c(z, divided.lnL[2])
    else
      z <- c(z, org.lnL[i])
  }
  z
}

# update the detVx by using the result of "split2cls()"
# continue: no update
# org.detVx: original detVx vector
# divided.detVx: divided detVx vector having 2 elements.
updtdetVx <- function(continue, org.detVx, k1, k2, divided.detVx){
  if (!is.vector(org.detVx))
    return(divided.detVx)
  if (!continue)
    return(org.detVx)
  if (k1 == k2)
    stop("updtdetVx() : k1 and k2 should differ.")
  
  z <- NULL
  for (i in 1:max(k2, length(org.detVx))){
    if (i == k1)
      z <- c(z, divided.detVx[1])
    else if (i == k2)
      z <- c(z, divided.detVx[2])
    else
      z <- c(z, org.detVx[i])
  }
  z
}

# split 2 clusters if we would prefer it based on BIC
# q: a number of parameters
# bic.prior: BIC which x is given; if bic.prior=NULL then we calculate
# lnL.prior: lnL which x is given; if bic.prior=NULL then we calculate
# detVx.prior: detVx which x is given; if bic.prior=NULL then we calculate
split2cls <- function(x, centers, q, bic.prior, lnL.prior, detVx.prior, iter.max, ignore.covar){
  if (is.null(bic.prior)){
    pb <- priorBIC(x, centers, q, ignore.covar)
    bic.prior <- pb$bic
    lnL.prior <- pb$lnL
    detVx.prior <- pb$detVx
  }
  bic.post <- postBICs(x, centers, q, iter.max, ignore.covar)
  
  subcluster <- bic.post$clsub$cluster
  #
  # compare whether if we should split
  if (is.na(bic.post$bic[3])){
    # BIC may has NA because of few data 
    continue <- FALSE
  }else if (bic.post$bic[3] < bic.prior){
    # splitting ...
    # replace the cluster number to cl$cluster
    continue <- TRUE
  }else{
    # not splitting...
    # return "subcluster" stored k1 
    continue <- FALSE
  }
  # note that "subcluster" gives 1 or 2 
  list(continue = continue, subcluster = subcluster, 
       bic.prior = bic.prior, bic.post = bic.post$bic,
       lnL.prior = lnL.prior, lnL.post = bic.post$lnL,
       detVx.prior = detVx.prior, detVx.post = bic.post$detVx,
       centers = bic.post$clsub$centers,
       clj1 = bic.post$clj1, clj2 = bic.post$clj2)
}




# return BIC (prior BIC)
priorBIC <- function(x, centers, q, ignore.covar){
  lnL0 <- lnL(x, centers, ignore.covar)
  bic <- -2 * lnL0$lnL + q * log(nrow(x)) # BIC
  # bic <- -2 * lnL0$lnL + q  # AIC
  list(lnL = lnL0$lnL, detVx = lnL0$detVx, bic = bic)
}


# return BICs (two posterior BICs)
postBICs <- function(x, centers, q, iter.max, ignore.covar){
  #
  # split to 2 clusters
  clsub <- kmeans(x, 2, iter.max)
  y.ok1 <- lapply(clsub$cluster, "==", 1) # 1st sub-cluster or not
  y.ok2 <- lapply(clsub$cluster, "==", 2) # 2nd sub-cluster or not
  # extract sub data
  p <- ncol(x)
  clj1 <- matrix(x[as.logical(y.ok1)], ncol=p)
  clj2 <- matrix(x[as.logical(y.ok2)], ncol=p)
  # ratio for pdf.
  r1 <- clsub$size[1] / sum(clsub$size)	# [0,1]
  r2 <- 1 - r1 	# [0,1]
  # two later BICs
  # print(clsub$centers[1,])	# for debug
  # print(apply(clj1,2,mean))	# for debug
  # print(sqrt(apply(clj1,2,var)))	# for debug
  # print(r1)	# for debug
  lnL1 <-  lnL(clj1, clsub$centers[1,], ignore.covar)
  # print(clsub$centers[2,])	# for debug
  # print(apply(clj2,2,mean))	# for debug
  # print(sqrt(apply(clj2,2,var)))	# for debug
  # print(r2)	# for debug
  lnL2 <-  lnL(clj2, clsub$centers[2,], ignore.covar)
  n1 <- nrow(clj1)
  n2 <- nrow(clj2)
  # normalizing factor; dist() is in library(mva)
  if (is.na(lnL1$detVx) || is.na(lnL2$detVx))
    beta <- 0
  else
    beta <- dist(clsub$center) / (sqrt(lnL1$detVx + lnL2$detVx))
  alpha <- 0.5 / pnorm(beta)
  BIC1 <- -2 * lnL1$lnL +q * log(n1)
  BIC2 <- -2 * lnL2$lnL +q * log(n2) 
  # BIC1 <- -2 * lnL1$lnL +q # AIC
  # BIC2 <- -2 * lnL2$lnL +q # AIC
  
  # cat (paste("alpha =",alpha,"\n"))	# for debug
  # cat (paste("beta =",beta,"\n"))	# for debug
  
  # BIC is not (BIC1 + BIC2)
  BIC <- -2 * lnL1$lnL  -2 * lnL2$lnL + 2 * q * log(n1 + n2) - 2 * (n1 + n2) * log(alpha)
  # BIC <- -2 * lnL1$lnL  -2 * lnL2$lnL + 2 * q  - 2 * (n1 + n2) * log(alpha) # AIC
  list(bic = c(BIC1, BIC2, BIC), 
       lnL = c(lnL1$lnL, lnL2$lnL),
       detVx = c(lnL1$detVx, lnL2$detVx),
       clsub = clsub, clj1 = clj1, clj2 = clj2)
}



# return BICs for Two-merged clusters model and devided clusters model
# k1/k2: marged cluster ID
mergedBIC <- function(x, xcl, k1, k2, q, ignore.covar, pr.proc){
  # sample size
  # check for input data
  n1 <- xcl$size[k1]
  n2 <- xcl$size[k2]
  if (n1 == 0 || n2 == 0){
    # already had been merged
    cat(paste("already had been merged\n"))
    ret <- F
    return( list (ret = ret))
  }
  if (is.null(xcl$lnL[k1]) || is.null(xcl$lnL[k2])){
    # lnL may be null because of few data
    cat(paste("lnL may be null because of few data\n"))
    ret <- F
    return( list (ret = ret))
  }
  
  # divided clusters model
  lnL1 = xcl$lnL[k1]
  lnL2 = xcl$lnL[k2]
  ctrextrt <- rbind(xcl$centers[k1,], xcl$centers[k2,])
  beta <- dist(ctrextrt) / (sqrt(xcl$detVx[k1] + xcl$detVx[k2]))
  if (pr.proc) cat(paste("beta=", round (beta, digit=2), "\n"))
  
  # if (beta > 10){
  # 	# 2 clusters far apart
  # 	ret <- F
  # 	return( list (ret = ret))
  # }
  
  alpha <- 0.5 / as.numeric(pnorm(beta))
  bicdiv <- -2 * lnL1  -2 * lnL2 + 2 * q * log(n1 + n2) - 2 * (n1 + n2) * log(alpha)
  # bicdiv <- -2 * lnL1 -2 * lnL2 + 2 * q - 2 * (n1 + n2) * log(alpha) # AIC
  
  # extract 2 clusters data
  y.ok1 <- lapply(xcl$cluster, "==", k1) # 1st sub-cluster or not
  y.ok2 <- lapply(xcl$cluster, "==", k2) # 2nd sub-cluster or not
  
  # extract sub data
  p = ncol(x)
  clj1 <- matrix(x[as.logical(y.ok1)], ncol=p)
  clj2 <- matrix(x[as.logical(y.ok2)], ncol=p)
  xmgd <- rbind(clj1, clj2)
  
  # merged cluster center
  ctrmgd <- (n1 * xcl$centers[k1,] + n2 * xcl$centers[k2,]) / (n1 + n2)
  lnLmgd <- lnL(xmgd, ctrmgd, ignore.covar)
  bicmgd <- -2 * lnLmgd$lnL + q * log(nrow(xmgd)) # BIC
  # bicmgd <- -2 * lnLmgd$lnL + q  # AIC
  
  ret <- T
  list (ret = ret, ctrmgd = ctrmgd, lnLmgd = lnLmgd$lnL, detVxmgd = lnLmgd$detVx, bicmgd = bicmgd, bicdiv = bicdiv)
}





# log-likelihood under the assumption of 
# p-dimensional multivariate normal distribution.
# ignore.covar: ignore the covariance 
lnL <- function(x, centers, ignore.covar=T){
  x <- as.matrix(x)
  p <- ncol(x)	# p-dimensional multivariate
  n <- nrow(x)	# sample size
  if (missing(centers)) 
    stop("centers must be a number or a matrix")
  if (n <= 2)	# few data
    return(list(lnL=NA, detVx=NA))
  vx <- var(x)	# var-co.var matrix
  # print(x)	# for debug
  if (p == 1){ # x is vector 
    invVx <- 1 / as.vector(vx)
    detVx <- as.vector(vx)
  }else{ 
    if (ignore.covar){
      invVx <- diag(1/diag(vx)) # inv. matrix when assuming diag.  
      detVx <- prod(diag(vx)) # det. when assuming diag. 
    }else{
      invVx <- solve(vx) # inverse matrix of "vx"
      y <- chol(vx) # Cholesky decomposition
      detVx <- prod(diag(y)) # vx = t(y) %*% y, where y is triangular,
      # then, det(vx) = det(t(y)) * det(y)
    }
  }
  t1 <- -p/2 * 1.837877066 # 1.837... = log(2 * 3.1415...)
  t2 <- -log(detVx) / 2
  xmu <- t(apply(x, 1, "-", centers))
  # print(centers)	# for debug
  # print(xmu)	# for debug
  # s <- 0
  # for (i in 1:n)
  #	s <- s + t(xmu[i,]) %*% invVx %*% xmu[i,]
  if (p == 1){
    s <- sum(xmu^2 * invVx)
  }else{
    s <- sum(apply(xmu, 1, txInvVxX, invVx=invVx))
  }
  t3 <- -s / 2
  ll <- (t1 + t2) * n + as.numeric(t3)	# log likelihood
  list(lnL=ll, detVx=detVx)
}

# function for calculation of 
# t(xmu[i,]) %*% invVx %*% xmu[i,]
txInvVxX <- function(x, invVx){
  t(x) %*% invVx %*% x
}