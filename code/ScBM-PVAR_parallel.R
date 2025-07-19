
get_sim <- function(q,TT,k,l,iter_max){

  s <- 4; p <- 1
  bd <- c(0.5,0.7,0.6,0.5)
  if (l == 1){
    bu <- 0.01; bl <- 0.01
  }else{
    bu <- 0.02; bl <- 0.02
  }
  
  B_list <- list()
  if (k == 1){
    nc <- c(2,2,2,2)
    n_comm <- c(2,2,2,2,2,2,2,2)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,nc[m]*nc[m]),nc[m],nc[m])
      diag(B_list[[m]]) <- rep(bd[m],nc[m])
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }else if (k == 2){
    nc <- c(4,4,4,4)
    n_comm <- c(4,4,4,4,4,4,4,4)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,nc[m]*nc[m]),nc[m],nc[m])
      diag(B_list[[m]]) <- rep(bd[m],nc[m])
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }else if (k == 3){
    nc <- c(2,3,3,2)
    n_comm <- c(2,3,3,3,3,3,3,2)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,n_comm[2*(m-1)+1]*n_comm[2*m]),n_comm[2*(m-1)+1],n_comm[2*m])
      if (length(diag(B_list[[m]])) < length(rep(bd[m],n_comm[2*(m-1)+1]))){
        diag(B_list[[m]]) <- rep(bd[m],n_comm[2*m])
      }else{
        diag(B_list[[m]]) <- rep(bd[m],n_comm[2*(m-1)+1])
      }
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }else {
    nc <- c(2,3,4,2)
    n_comm <- c(2,3,3,4,4,4,4,2)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,n_comm[2*(m-1)+1]*n_comm[2*m]),n_comm[2*(m-1)+1],n_comm[2*m])
      if (length(diag(B_list[[m]])) < length(rep(bd[m],n_comm[2*(m-1)+1]))){
        diag(B_list[[m]]) <- rep(bd[m],n_comm[2*m])
      }else{
        diag(B_list[[m]]) <- rep(bd[m],n_comm[2*(m-1)+1])
      }
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }
  
  library(doParallel)
  require(dqrng)
  require(doRNG)

  cl <- makeCluster(detectCores())
  registerDoParallel(cl)

  PVAR_sim <- foreach(iter = 1:iter_max,
                      .errorhandling = 'remove',
                      .packages = c("igraph","MASS","cluster","combinat","foreach",
                                    "mclust","gtools")) %dopar% {

                        source("ScBM_library.R")
                        s <- 4; p <- 1
                        Sigma <- diag(0.1,q)
                        
                        # --------------------------------------------------------------- #
                        # Data generation:
                        phi <- list()
                        for (m in 1:s){
                          phi[[m]] <- 0.99
                        }
                        
                        max_eig <- 1
                        while (max_eig > 0.95){
                          ScBM <- Phi_generator(q, B_list, s, phi, selfTF=TRUE, 
                                                n_comm, threshold=TRUE)
                          
                          Phi_large <- matrix(0,s*q,s*q)
                          Phi_large[1:q,(3*q+1):(4*q)] <- ScBM$Phi[[1]]
                          Phi_large[(q+1):(2*q),1:q] <- ScBM$Phi[[2]]
                          Phi_large[(2*q+1):(3*q),(q+1):(2*q)] <- ScBM$Phi[[3]]
                          Phi_large[(3*q+1):(4*q),(2*q+1):(3*q)] <- ScBM$Phi[[4]]
                          
                          max_eig <- spectral_radius_power(Phi_large)
                          if (max_eig > 0.95){
                            cat("Again!\n")
                            phi[[1]] <- phi[[1]]*0.9
                            phi[[2]] <- phi[[2]]*0.9
                            phi[[3]] <- phi[[3]]*0.9
                            phi[[4]] <- phi[[4]]*0.9
                          }
                        }
                        Phi_series <- cbind(ScBM$Phi[[1]],ScBM$Phi[[2]],ScBM$Phi[[3]],ScBM$Phi[[4]])
                        Yt <- PVAR_generator(TT, s, p, Phi_series, Sigma=diag(1,q))
                        Phi_hat <- PVAR_ols(Yt,s,p=1)$Phi_hat

                        # --------------------------------------------------------------- #
                        # alpha + kmeans
                        alpha_hat_kmeans <- alpha_cv(Phi_hat,q,p,n_comm,type="PVAR",clst="kmeans",fold=5)

                        PisCES_hat_kmeans <- PisCES(Phi_hat,q,p,n_comm,s,alpha=alpha_hat_kmeans$alpha_hat)
                        PVAR_4R1L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[4]],
                                                       PisCES_hat_kmeans$VL_bar[[1]]),q,nc[1],clst="kmeans")
                        PVAR_1R2L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[1]],
                                                       PisCES_hat_kmeans$VL_bar[[2]]),q,nc[2],clst="kmeans")
                        PVAR_2R3L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[2]],
                                                       PisCES_hat_kmeans$VL_bar[[3]]),q,nc[3],clst="kmeans")
                        PVAR_3R4L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[3]],
                                                       PisCES_hat_kmeans$VL_bar[[4]]),q,nc[4],clst="kmeans")

                        output_kmeans <- list(evaluate_clustering(ScBM$z_vec[,1],PVAR_1R2L_kmeans$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,2],PVAR_2R3L_kmeans$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,3],PVAR_3R4L_kmeans$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,4],PVAR_4R1L_kmeans$group_1q))

                        # --------------------------------------------------------------- #
                        # alpha + pam
                        alpha_hat_pam <- alpha_cv(Phi_hat,q,p,n_comm,type="PVAR",clst="pam",fold=5)

                        PisCES_hat_pam <- PisCES(Phi_hat,q,p,n_comm,s,alpha=alpha_hat_pam$alpha_hat)
                        PVAR_4R1L_pam <- Blockbuster(cbind(PisCES_hat_pam$VR_bar[[4]],
                                                       PisCES_hat_pam$VL_bar[[1]]),q,nc[1],clst="pam")
                        PVAR_1R2L_pam <- Blockbuster(cbind(PisCES_hat_pam$VR_bar[[1]],
                                                       PisCES_hat_pam$VL_bar[[2]]),q,nc[2],clst="pam")
                        PVAR_2R3L_pam <- Blockbuster(cbind(PisCES_hat_pam$VR_bar[[2]],
                                                       PisCES_hat_pam$VL_bar[[3]]),q,nc[3],clst="pam")
                        PVAR_3R4L_pam <- Blockbuster(cbind(PisCES_hat_pam$VR_bar[[3]],
                                                       PisCES_hat_pam$VL_bar[[4]]),q,nc[4],clst="pam")

                        output_pam <- list(evaluate_clustering(ScBM$z_vec[,1],PVAR_1R2L_pam$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,2],PVAR_2R3L_pam$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,3],PVAR_3R4L_pam$group_1q),
                                              evaluate_clustering(ScBM$z_vec[,4],PVAR_4R1L_pam$group_1q))

                        # --------------------------------------------------------------- #
                        # no alpha + kmeans
                        PisCES_hat_noalpha <- PisCES(Phi_hat,q,p,n_comm,s,alpha=0)
                        PVAR_4R1L_noalpha <- Blockbuster(cbind(PisCES_hat_noalpha$VR_bar[[4]],
                                                        PisCES_hat_noalpha$VL_bar[[1]]),q,nc[1],clst="kmeans")
                        PVAR_1R2L_noalpha <- Blockbuster(cbind(PisCES_hat_noalpha$VR_bar[[1]],
                                                        PisCES_hat_noalpha$VL_bar[[2]]),q,nc[2],clst="kmeans")
                        PVAR_2R3L_noalpha <- Blockbuster(cbind(PisCES_hat_noalpha$VR_bar[[2]],
                                                        PisCES_hat_noalpha$VL_bar[[3]]),q,nc[3],clst="kmeans")
                        PVAR_3R4L_noalpha <- Blockbuster(cbind(PisCES_hat_noalpha$VR_bar[[3]],
                                                        PisCES_hat_noalpha$VL_bar[[4]]),q,nc[4],clst="kmeans")
                        output_noalpha <- list(evaluate_clustering(ScBM$z_vec[,1],PVAR_1R2L_noalpha$group_1q),
                                           evaluate_clustering(ScBM$z_vec[,2],PVAR_2R3L_noalpha$group_1q),
                                           evaluate_clustering(ScBM$z_vec[,3],PVAR_3R4L_noalpha$group_1q),
                                           evaluate_clustering(ScBM$z_vec[,4],PVAR_4R1L_noalpha$group_1q))

                        # --------------------------------------------------------------- #
                        # no alpha + pam
                        PisCES_hat_noalpha_pam <- PisCES(Phi_hat,q,p,n_comm,s,alpha=0)
                        PVAR_4R1L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[4]],
                                                         PisCES_hat_noalpha_pam$VL_bar[[1]]),q,nc[1],clst="pam")
                        PVAR_1R2L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[1]],
                                                         PisCES_hat_noalpha_pam$VL_bar[[2]]),q,nc[2],clst="pam")
                        PVAR_2R3L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[2]],
                                                         PisCES_hat_noalpha_pam$VL_bar[[3]]),q,nc[3],clst="pam")
                        PVAR_3R4L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[3]],
                                                         PisCES_hat_noalpha_pam$VL_bar[[4]]),q,nc[4],clst="pam")
                        output_noalpha_pam <- list(evaluate_clustering(ScBM$z_vec[,1],PVAR_1R2L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,2],PVAR_2R3L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,3],PVAR_3R4L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,4],PVAR_4R1L_noalpha_pam$group_1q))
                        

                        # --------------------------------------------------------------- #
                        # output
                        output <- list()
                        output$Phi_hat <- Phi_hat
                        output$y_vec <- ScBM$y_vec
                        output$z_vec <- ScBM$z_vec

                        output$kmeans <- output_kmeans
                        output$pam <- output_pam
                        output$noalpha <- output_noalpha
                        output$noalpha_pam <- output_noalpha_pam
                        output
  }
  stopCluster(cl)
  print(length(PVAR_sim))
  if (length(PVAR_sim) == 0) {
    warning("No successful iterations were returned.")
  }
  save(PVAR_sim,file=paste0("sim_PVAR_q",q,"_T",TT,"_B_mat_type",k,"_deg_type",l,".RData"))
}

#####################################################################################

# Run iterations
q_list <- c(24,60,120)
T_list <- c(200,400,1000,2000,4000)
for (i in 1:3){
  q <- q_list[i]

  for (j in 1:5){
    TT <- T_list[j]

    for (k in 1:4){

      for (l in 1:2){

        cat("The current iteration is q:",q,"TT:",TT,"B_mat_type:",k,"deg_type:",l,"\n")
        get_sim(q,TT,k,l,iter_max=100)
        cat("Done.\n")
      }
    }
  }
}

######################################################################################

# Summarize the results
library(writexl)
q_list <- c(24,60,120)
T_list <- c(200,400,1000,2000,4000)

PVAR_summary <- data.frame()
cnt <- 1
for (i in 1:3){
  q <- q_list[i]

  for (j in 1:5){
    TT <- T_list[j]

    for (k in 1:4){

      for (l in 1:2){

        cat("The current iteration is q:",q,"TT:",TT,"B_mat_type:",k,"deg_type:",l,"\n")
        load(paste0("sim_PVAR_q",q,"_T",TT,"_B_mat_type",k,"_deg_type",l,".RData"))

        iter_num <- length(PVAR_sim)

        # kmeans:
        acc_kmeans <- round(1 - mean(mapply(iter=1:iter_num,
                              function(iter)mean(mapply(m=1:4,
                              function(m)PVAR_sim[[iter]]$kmeans[[m]]$err[1])))),3)
        ARI_kmeans <- round(mean(mapply(iter=1:iter_num,
                                  function(iter)mean(mapply(m=1:4,
                                  function(m)PVAR_sim[[iter]]$kmeans[[m]]$err[2])))),3)

        # pam:
        acc_pam <- round(1 - mean(mapply(iter=1:iter_num,
                                      function(iter)mean(mapply(m=1:4,
                                      function(m)PVAR_sim[[iter]]$pam[[m]]$err[1])))),3)
        ARI_pam <- round(mean(mapply(iter=1:iter_num,
                                  function(iter)mean(mapply(m=1:4,
                                  function(m)PVAR_sim[[iter]]$pam[[m]]$err[2])))),3)

        # no alpha:
        acc_noalpha <- round(1 - mean(mapply(iter=1:iter_num,
                                   function(iter)mean(mapply(m=1:4,
                                   function(m)PVAR_sim[[iter]]$noalpha[[m]]$err[1])))),3)
        ARI_noalpha <- round(mean(mapply(iter=1:iter_num,
                               function(iter)mean(mapply(m=1:4,
                               function(m)PVAR_sim[[iter]]$noalpha[[m]]$err[2])))),3)

        # no alpha + pam:
        acc_noalpha_pam <- round(1 - mean(mapply(iter=1:iter_num,
                                   function(iter)mean(mapply(m=1:4,
                                   function(m)PVAR_sim[[iter]]$noalpha_pam[[m]]$err[1])))),3)
        ARI_noalpha_pam <- round(mean(mapply(iter=1:iter_num,
                                   function(iter)mean(mapply(m=1:4,
                                   function(m)PVAR_sim[[iter]]$noalpha_pam[[m]]$err[2])))),3)

        PVAR_summary <- rbind(PVAR_summary,data.frame(q=q,
                                                      T=TT,
                                                      scn=k,
                                                      off=l,
                                                      acc_kmeans=acc_kmeans,
                                                      acc_noalpha=acc_noalpha,
                                                      acc_pam=acc_pam,
                                                      acc_noalpha_pam=acc_noalpha_pam,
                                                      ARI_kmeans=ARI_kmeans,
                                                      ARI_noalpha=ARI_noalpha,
                                                      ARI_pam=ARI_pam,
                                                      ARI_noalpha_pam=ARI_noalpha_pam))

        cnt <- cnt + 1
        cat("Done.\n")
      }
    }
  }
}
write_xlsx(PVAR_summary, "PVAR_summary.xlsx")


