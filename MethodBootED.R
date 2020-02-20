# Bootstrap method, size of bootstrap = 1
# Update Dec 7, 2019


MethodBoot <- function(tree, form1, n_boot, boot_r=1){
  data <- tree$data
  PID <- c(1:nrow(data))
  data <-cbind(data,PID)
  ncol<-ncol(data)
  model_formula <- form1
  nleaves_orig <- nrow(tree$li)
  # Person ID number and Condition (which is always in the second column of the data set)
  PID_Cond_orig <- as.data.frame(cbind(data[, 2],PID))
  leaves_origTree_G1 <- list()
  leaves_origTree_G2 <- list()

  for (i in (1:nleaves_orig)){
    leaf_index <- tree$nind[,i]
    leaf_assign_orig <- subset(cbind(PID_Cond_orig, leaf_index), leaf_index == 1)
    LiG1 <- subset( leaf_assign_orig,leaf_assign_orig[,1] == '1')
    LiG2 <- subset(leaf_assign_orig, leaf_assign_orig[,1] == '2')
    leaves_origTree_G1[[i]] <- LiG1
    leaves_origTree_G2[[i]] <- LiG2
  }

  #print(paste('Subset 1', leaves_origTree_G1[1]))
  boot_matrix_G1 <- matrix(nrow = n_boot, ncol = nleaves_orig)
  boot_matrix_G2 <- matrix(nrow = n_boot, ncol = nleaves_orig)
  boot_matrix_d <- matrix(nrow = n_boot, ncol = nleaves_orig)
  boot_nleaves <- numeric(n_boot)
 #Start bootstrapping
  for (n in 1: n_boot){
    bsample <- data[sample(nrow(data), boot_r*nrow(data), replace = TRUE), ]
    #print(head(sample))
    control2 <- tree$control
    T_samp_unpruned <- quint(model_formula, bsample, control=control2)
    T_samp <- prune(T_samp_unpruned)
    #plot(T_samp)
    nleaves_samp <- nrow(T_samp$li)
    boot_nleaves[n] <- nleaves_samp #save the number of leaves of the pruned tree in the bootstrap sample
    PID_Cond_samp <- as.data.frame(bsample[, c(2,ncol)])
    #create a cross-table: rows are the number of leaves of the bootstrap tree
    #and columns are number of leaves in original tree
     matrix_treat1 <- matrix(nrow = nleaves_samp, ncol = nleaves_orig)
    matrix_treat2 <- matrix(nrow = nleaves_samp, ncol = nleaves_orig)
    matrixn_treat1 <- matrix(nrow = nleaves_samp, ncol = nleaves_orig)
    matrixn_treat2 <- matrix(nrow = nleaves_samp, ncol = nleaves_orig)

    for (i in 1:nleaves_samp){       #nleaves_samp
      #print(i)
      for (j in 1:nleaves_orig){     #nleaves_orig
        #print(j)
        leaf_index <- T_samp$nind[,i]
        leaf_assign_sample <- subset(cbind(PID_Cond_samp, leaf_index), leaf_index == 1)

        LiG1_samp <- subset(leaf_assign_sample,leaf_assign_sample[,1]  == '1')
        LiG2_samp <- subset(leaf_assign_sample, leaf_assign_sample[,1] == '2')
        LjG1_orig <- leaves_origTree_G1[[j]]
        LjG2_orig <- leaves_origTree_G2[[j]]

        nzij_1 <- length(LiG1_samp$PID[is.element(LiG1_samp$PID,LjG1_orig$PID)])
        nzij_2 <- length(LiG2_samp$PID[is.element(LiG2_samp$PID,LjG2_orig$PID)])
       #nzij_1: the intersection in cell ij = number of the persons in cell ij in treatment 1

        mean_bootij_1 <- T_samp$li$`meanY|T=1`[i] #bootstrap mean of treatment 1 in Leaf 1
        mean_bootij_2 <- T_samp$li$`meanY|T=2`[i] #bootstrap mean of treatment 2 in Leaf 1
        #print(paste('Mean_1:', mean_bootij_1))
        cell_productij_1 <- nzij_1 * mean_bootij_1
        cell_productij_2 <- nzij_2 * mean_bootij_2

        matrixn_treat1[i,j] <- nzij_1
        matrixn_treat2[i,j] <- nzij_2
        matrix_treat1[i,j] <- cell_productij_1
        matrix_treat2[i,j] <- cell_productij_2
      }
    }
    # matrixn_treat1: the rows add up to the number of persons in treatment 1 in the bootstrap tree
    mu_G1 <- numeric(nleaves_orig)
    mu_G2 <- numeric(nleaves_orig)
    d_mu <- numeric(nleaves_orig)
    for (i in 1:ncol(matrix_treat1)){
      mu_G1[i] <- sum(matrix_treat1[,i]) / sum(matrixn_treat1[,i])
      mu_G2[i]<- sum(matrix_treat2[,i]) / sum(matrixn_treat2[,i])
      d_mu[i] <- mu_G1[i] -  mu_G2[i]
  }
    boot_matrix_G1[n,] <- mu_G1
    boot_matrix_G2[n,] <- mu_G2
    boot_matrix_d[n, ] <- d_mu
  }
  obj <- list(nleaves= boot_nleaves, meanT_1 = boot_matrix_G1, meanT_2 = boot_matrix_G2,
              meandif=boot_matrix_d)
  return(obj)
}
