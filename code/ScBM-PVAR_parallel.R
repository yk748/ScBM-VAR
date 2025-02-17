
get_sim <- function(q,TT,k,l,iter_max){
  
    if (k == 1){
      nc <- 2
    }else if (k == 2){
      nc <- 3
    }else if (k == 3){
      nc <- 5
    }

    if (l == 1){
      bd <- 0.5; bu <- 0.05; bl <- 0.1
    }else if (l == 2){
      bd <- 0.5; bu <- 0.15; bl <- 0.1
    }else if (l == 3){
      bd <- 0.5; bu <- 0.25; bl <- 0.1
    }
    
    K <- 5; s <- 4; p <- 1
    B_list <- list()
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,K*K),K,K)
      diag(B_list[[m]]) <- rep(bd,K)
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
      if (nc != K){
        B_list[[m]][,(nc+1):K] <- B_list[[m]][(nc+1):K,] <- 0
      }
    }
    
    library(doParallel)
    require(dqrng)
    require(doRNG)
    
    list_dummy <- list()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    
    PVAR_sim <- foreach(iter = 1:iter_max,
                        .errorhandling = 'remove',
                        .packages = c("MASS","combinat","foreach","randnet","mclust")) %dopar% {
                          
      source('ScBM-VAR_library.R')   
      K <- 5; s <- 4; p <- 1          
      Sigma <- diag(1,q)
      
      phi <- list() # Some specific constant
      phi[[1]] <- 0.99; phi[[2]] <- 0.99; phi[[3]] <- 0.99; phi[[4]] <- 0.99
      
      max_eig <- 1
      while (max_eig > 0.95){
        ScBM <- DC_ScBM_generator(q,s,rep(K,2*s),B_list,
                                  type="PVAR",p,phi=phi,tau=TRUE)
        Phi_large <- matrix(0,s*q,s*q)
        Phi_large[1:q,(3*q+1):(4*q)] <- ScBM$Phi[[1]]
        Phi_large[(q+1):(2*q),1:q] <- ScBM$Phi[[2]]
        Phi_large[(2*q+1):(3*q),(q+1):(2*q)] <- ScBM$Phi[[3]]
        Phi_large[(3*q+1):(4*q),(2*q+1):(3*q)] <- ScBM$Phi[[4]]
        
        max_eig <- max(abs(eigen(Phi_large)$values))
        if (max_eig > 0.95){
          cat("Again!\n")
          phi[[1]] <- phi[[1]]*0.9
          phi[[2]] <- phi[[2]]*0.9
          phi[[3]] <- phi[[3]]*0.9
          phi[[4]] <- phi[[4]]*0.9
        }
      }
      Phi_series <- cbind(ScBM$Phi[[1]],ScBM$Phi[[2]],ScBM$Phi[[3]],ScBM$Phi[[4]])
      Yt <- PVAR_generator(TT, s, p, Phi_series, Sigma)
      Phi_hat <- PVAR_ols_new(Yt,s,p=1)$Phi_hat
      
      alpha_hat <- alpha_cv(Phi_hat,q,p,rep(K,2*s),type="PVAR",fold=5)
      
      PisCES_hat <- PisCES(Phi_hat,q,p,rep(K,2*s),s,alpha=alpha_hat$alpha_hat)
      PVAR_4R1L <- Blockbuster(cbind(PisCES_hat$VR_bar[[4]],
                                     PisCES_hat$VL_bar[[1]]),q,K)
      PVAR_1R2L <- Blockbuster(cbind(PisCES_hat$VR_bar[[1]],
                                     PisCES_hat$VL_bar[[2]]),q,K)
      PVAR_2R3L <- Blockbuster(cbind(PisCES_hat$VR_bar[[2]],
                                     PisCES_hat$VL_bar[[3]]),q,K)
      PVAR_3R4L <- Blockbuster(cbind(PisCES_hat$VR_bar[[3]],
                                     PisCES_hat$VL_bar[[4]]),q,K)
      
      acc_agg <- array(NA,dim=c(s,2)); acc_class <- array(NA,dim=c(s,2))
      ARI_tab <- array(NA,dim=c(s,2)); NMI_tab <- array(NA,dim=c(s,2))
      for (m in 1:s){
        
        ref_L <- mapply(i=1:q,function(i) which(ScBM$Y[[m]][i,]==1))
        ref_R <- mapply(i=1:q,function(i) which(ScBM$Z[[m]][i,]==1))
        
        if (m == 1){
          est_L <- PVAR_4R1L$group_1q
          est_R <- PVAR_1R2L$group_1q
        }else if (m == 2){
          est_L <- PVAR_1R2L$group_1q
          est_R <- PVAR_2R3L$group_1q
        }else if (m == 3){
          est_L <- PVAR_2R3L$group_1q
          est_R <- PVAR_3R4L$group_1q
        }else if (m == 4){
          est_L <- PVAR_3R4L$group_1q
          est_R <- PVAR_4R1L$group_1q
        }
        perm_L <- lapply(permn(c(1:K)), function(purt) {
          new_L <- est_L
          for (i in 1:K) {
            new_L[est_L == i] <- purt[i]
          }
          new_L
        })
        perm_R <- lapply(permn(c(1:K)), function(purt) {
          new_R <- est_R
          for (i in 1:K) {
            new_R[est_R == i] <- purt[i]
          }
          new_R
        })
        
        acc_agg_L <- vector("numeric",factorial(K)) 
        acc_agg_R <- vector("numeric",factorial(K))
        acc_class_L <- vector("numeric",factorial(K)) 
        acc_class_R <- vector("numeric",factorial(K))
        for (i in 1:factorial(K)){
          acc_agg_L[i] <- sum(perm_L[[i]] == ref_L)/q
          acc_agg_R[i] <- sum(perm_R[[i]] == ref_R)/q
          acc_class_L[i] <- mean(colMeans(mapply(k=1:K,function(k)perm_L[[i]] == k)[,1:nc]
                                          == mapply(k=1:K,function(k)ref_L==k)[,1:nc]))
          acc_class_R[i] <- mean(colMeans(mapply(k=1:K,function(k)perm_R[[i]] == k)[,1:nc]
                                          == mapply(k=1:K,function(k)ref_R==k)[,1:nc]))
        }
        acc_agg[m,1] <- max(acc_agg_L); acc_agg[m,2] <- max(acc_agg_R);
        acc_class[m,1] <- max(acc_class_L); acc_class[m,2] <- max(acc_class_R);
        ARI_tab[m,1] <- adjustedRandIndex(perm_L[[1]],ref_L) ; 
        ARI_tab[m,2] <- adjustedRandIndex(perm_R[[1]],ref_R);
        NMI_tab[m,1] <- NMI(perm_L[[1]],ref_L); 
        NMI_tab[m,2] <- NMI(perm_R[[1]],ref_R);
      }
      
      output <- list()
      output$Y <- ScBM$Y
      output$Z <- ScBM$Z
      output$PisCES <- PisCES_hat
      output$PVAR_4R1L <- PVAR_4R1L
      output$PVAR_1R2L <- PVAR_1R2L
      output$PVAR_2R3L <- PVAR_2R3L
      output$PVAR_3R4L <- PVAR_3R4L
      output$acc_agg_L <- acc_agg_L
      output$acc_agg_R <- acc_agg_R
      output$acc_agg <- acc_agg
      output$acc_class_L <- acc_class_L
      output$acc_class_R <- acc_class_R
      output$acc_class <- acc_class
      output$ARI <- ARI_tab
      output$NMI <- NMI_tab
      list_dummy[[iter]] <- output
    }  
    stopCluster(cl)
    save(PVAR_sim,file=paste0("PVAR_q",q,"_T",TT,"_nc",nc,"_Bmat_type",l,".RData"))
}

#####################################################################################

q_list <- c(20,30,50,100)
T_list <- c(100,200,400,1000,2000,4000)
for (i in 1:4){
  q <- q_list[i]
  
  for (j in 1:6){
    TT <- T_list[j]

    for (k in c(1,3)){
      
      for (l in 1:3){
        cat("The current iteration is q:",q,"TT:",TT,"k:",k,"l:",l,"\n")
        get_sim(q,TT,k,l,iter_max=100)
        cat("Done.\n")
      }
    }
  }
}

#####################################################################################

q_list <- c(20,30,50,100)
T_list <- c(100,200,400,1000,2000,4000)

list_total <- list()
cnt <- 1
for (i in 1:4){
  q <- q_list[i]

  for (j in 1:6){
    TT <- T_list[j]

    for (k in 1:3){
      if (k == 1){
        nc <- 2
      }else if (k == 2){
        nc <- 3
      }else if (k == 3){
        nc <- 5
      }

      for (l in 1:3){

        cat("Current q is",q,"T is",TT,"nc is",nc,"Bmat is",l,"\n")
        load(paste0("PVAR_q",q,"_T",TT,"_nc",nc,"_Bmat_type",l,".RData"))

        iter_num <- length(PVAR_sim)
        tab_agg <- colMeans(matrix(rowMeans(mapply(x=1:iter_num,
                                          function(x)PVAR_sim[[x]]$acc_agg)),ncol=2))
        tab_class <- colMeans(matrix(rowMeans(mapply(x=1:iter_num,
                                            function(x)PVAR_sim[[x]]$acc_class)),ncol=2))
        tab_ARI <- colMeans(matrix(rowMeans(mapply(x=1:iter_num,
                                          function(x)PVAR_sim[[x]]$ARI)),ncol=2))
        tab_NMI <- colMeans(matrix(rowMeans(mapply(x=1:iter_num,
                                          function(x)PVAR_sim[[x]]$NMI)),ncol=2))

        list_total[[cnt]] <- list(q=q,TT=TT,nc=nc,Bmat=l,counts=iter_num,
                                  agg=tab_agg,class=tab_class,ARI=tab_ARI,NMI=tab_NMI)
        cnt <- cnt + 1
      }
    }
  }
}
save(list_total,file="PVAR_sim_total.RData")


#####################################################################################
load("PVAR_sim_total.RData")

for (i in 1:216){
  cat("#--------------------------#\n")
  cat("Current q is",list_total[[i]]$q,"T is",list_total[[i]]$TT,"nc is",list_total[[i]]$nc,"and Bmat is",list_total[[i]]$Bmat,"\n")
  cat(c(round(list_total[[i]]$class[1],3),"&",round(list_total[[i]]$NMI[1],3)),"\n")
}


