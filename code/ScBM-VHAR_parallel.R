
get_sim <- function(q,TT,k,l,iter_max){

  s <- 3; p <- 1
  bd <- c(0.7,0.8,0.9)
  if (l == 1){
    bu <- 0.01; bl <- 0.01
  }else{
    bu <- 0.02; bl <- 0.02
  }
  
  B_list <- list()
  if (k == 1){
    nc <- c(2,2,2)
    n_comm <- c(2,2,2,2,2,2)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,nc[m]*nc[m]),nc[m],nc[m])
      diag(B_list[[m]]) <- rep(bd[m],nc[m])
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }else if (k == 2){
    nc <- c(4,4,4)
    n_comm <- c(4,4,4,4,4,4)
    
    for (m in 1:s){
      B_list[[m]] <- matrix(rep(0,nc[m]*nc[m]),nc[m],nc[m])
      diag(B_list[[m]]) <- rep(bd[m],nc[m])
      B_list[[m]][upper.tri(B_list[[m]])] <- bu
      B_list[[m]][lower.tri(B_list[[m]])] <- bl
    }
    
  }else if (k == 3){
    nc <- c(2,2,3)
    n_comm <- c(2,2,2,3,3,3)
    
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
    nc <- c(4,2,2)
    n_comm <- c(4,4,4,2,2,2)
    
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

  VHAR_sim <- foreach(iter = 1:iter_max,
                      .errorhandling = 'remove',
                      .packages = c("igraph","MASS","cluster","combinat","foreach",
                                    "mclust","gtools")) %dopar% {

                        source("ScBM_library.R")
                        s <- 3; p <- 1
                        Sigma <- diag(1,q)

                        # --------------------------------------------------------------- #
                        # Data generation:
                        phi <- list()
                        for (m in 1:s){
                          phi[[m]] <- 0.8
                        }

                        max_eig <- 1
                        while(max_eig > 0.95){
                          ScBM <- Phi_generator(q, B_list, s, phi, selfTF=TRUE,n_comm,threshold=TRUE)
                          Phi_mat <- array(NA,dim=c(q,q,22))
                          Phi_mat[,,1] <- ScBM$Phi[[1]]+1/5*ScBM$Phi[[2]]+1/22*ScBM$Phi[[3]]
                          Phi_mat[,,2] <- Phi_mat[,,3] <- Phi_mat[,,4] <- Phi_mat[,,5] <- 1/5*ScBM$Phi[[2]]+1/22*ScBM$Phi[[3]]
                          for (h in 6:22){
                            Phi_mat[,,h] <- 1/22*ScBM$Phi[[3]]
                          }
                          Phi_large <- companion_form_phi(Phi_mat,d=q,p=22)

                          max_eig <- spectral_radius_power(Phi_large)
                          if (max_eig > 0.95){
                            cat("Again!\n")
                            phi[[1]] <- phi[[1]]*0.9
                            phi[[2]] <- phi[[2]]*0.9
                            phi[[3]] <- phi[[3]]*0.9
                          }
                        }
                        Phi_series <- cbind(ScBM$Phi[[1]],ScBM$Phi[[2]],ScBM$Phi[[3]])
                        Yt <- VHAR_generator(TT,Phi_series,Sigma=diag(1,q))
                        Phi_hat <- VHAR_ols(Yt)$Phi_hat

                        # --------------------------------------------------------------- #
                        # alpha + kmeans
                        alpha_hat_kmeans <- alpha_cv(Phi_hat,q,p,n_comm,type="VHAR",clst="kmeans",fold=5)

                        PisCES_hat_kmeans <- PisCES(Phi_hat,q,p,n_comm,s,alpha=alpha_hat_kmeans$alpha_hat)
                        VHAR_1L_kmeans <- Blockbuster(PisCES_hat_kmeans$VL_bar[[1]],q,nc[1],clst="kmeans")
                        VHAR_1R2L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[1]],
                                                       PisCES_hat_kmeans$VL_bar[[2]]),q,nc[2],clst="kmeans")
                        VHAR_2R3L_kmeans <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[2]],
                                                       PisCES_hat_kmeans$VL_bar[[3]]),q,nc[2],clst="kmeans")
                        VHAR_3R_kmeans <- Blockbuster(PisCES_hat_kmeans$VR_bar[[3]],q,nc[3],clst="kmeans")
                        
                        output_kmeans <- list(evaluate_clustering(ScBM$y_vec[,1],VHAR_1L_kmeans$group_1q),
                                           evaluate_clustering(ScBM$y_vec[,2],VHAR_1R2L_kmeans$group_1q),
                                           evaluate_clustering(ScBM$y_vec[,3],VHAR_2R3L_kmeans$group_1q),
                                           evaluate_clustering(ScBM$z_vec[,3],VHAR_3R_kmeans$group_1q))

                        # --------------------------------------------------------------- #
                        # alpha + pam
                        alpha_hat_pam <- alpha_cv(Phi_hat,q,p,n_comm,type="VHAR",clst="pam",fold=5)

                        PisCES_hat_pam <- PisCES(Phi_hat,q,p,n_comm,s,alpha=alpha_hat_pam$alpha_hat)
                        VHAR_1L_pam <- Blockbuster(PisCES_hat_kmeans$VL_bar[[1]],q,nc[1],clst="pam")
                        VHAR_1R2L_pam <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[1]],
                                                          PisCES_hat_kmeans$VL_bar[[2]]),q,nc[2],clst="pam")
                        VHAR_2R3L_pam <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[2]],
                                                          PisCES_hat_kmeans$VL_bar[[3]]),q,nc[2],clst="pam")
                        VHAR_3R_pam <- Blockbuster(PisCES_hat_kmeans$VR_bar[[3]],q,nc[3],clst="pam")
                        
                        output_pam <- list(evaluate_clustering(ScBM$y_vec[,1],VHAR_1L_pam$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,2],VHAR_1R2L_pam$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,3],VHAR_2R3L_pam$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,3],VHAR_3R_pam$group_1q))

                        # --------------------------------------------------------------- #
                        # no alpha + kmeans
                        PisCES_hat_noalpha <- PisCES(Phi_hat,q,p,rep(nc,2*s),s,alpha=0)
                        VHAR_1L_noalpha <- Blockbuster(PisCES_hat_kmeans$VL_bar[[1]],q,nc[1],clst="kmeans")
                        VHAR_1R2L_noalpha <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[1]],
                                                       PisCES_hat_kmeans$VL_bar[[2]]),q,nc[2],clst="kmeans")
                        VHAR_2R3L_noalpha <- Blockbuster(cbind(PisCES_hat_kmeans$VR_bar[[2]],
                                                       PisCES_hat_kmeans$VL_bar[[3]]),q,nc[2],clst="kmeans")
                        VHAR_3R_noalpha <- Blockbuster(PisCES_hat_kmeans$VR_bar[[3]],q,nc[3],clst="kmeans")

                        
                        output_noalpha <- list(evaluate_clustering(ScBM$y_vec[,1],VHAR_1L_noalpha$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,2],VHAR_1R2L_noalpha$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,3],VHAR_2R3L_noalpha$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,3],VHAR_3R_noalpha$group_1q))
                        
                        # --------------------------------------------------------------- #
                        # no alpha + pam
                        PisCES_hat_noalpha_pam <- PisCES(Phi_hat,q,p,rep(nc,2*s),s,alpha=0)
                        VHAR_1L_noalpha_pam <- Blockbuster(PisCES_hat_noalpha_pam$VL_bar[[1]],q,nc[1],clst="pam")
                        VHAR_1R2L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[1]],
                                                               PisCES_hat_noalpha_pam$VL_bar[[2]]),q,nc[2],clst="pam")
                        VHAR_2R3L_noalpha_pam <- Blockbuster(cbind(PisCES_hat_noalpha_pam$VR_bar[[2]],
                                                               PisCES_hat_noalpha_pam$VL_bar[[3]]),q,nc[2],clst="pam")
                        VHAR_3R_noalpha_pam <- Blockbuster(PisCES_hat_noalpha_pam$VR_bar[[3]],q,nc[3],clst="pam")
                        
                        
                        output_noalpha_pam <- list(evaluate_clustering(ScBM$y_vec[,1],VHAR_1L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,2],VHAR_1R2L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$y_vec[,3],VHAR_2R3L_noalpha_pam$group_1q),
                                               evaluate_clustering(ScBM$z_vec[,3],VHAR_3R_noalpha_pam$group_1q))

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
  print(length(VHAR_sim))
  if (length(VHAR_sim) == 0) {
    warning("No successful iterations were returned.")
  }
  save(VHAR_sim,file=paste0("sim_VHAR_q",q,"_T",TT,"_B_mat_type",k,"_deg_type",l,".RData"))
}

#####################################################################################

# Run iterations
q_list <- c(24,60,120)
T_list <- c(400,1000,2000,4000,8000)
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

# ######################################################################################

# Summarize the results
library(writexl)
q_list <- c(24,60,120)
T_list <- c(400,1000,2000,4000,8000)

VHAR_summary <- data.frame()
cnt <- 1
for (i in 1:3){
  q <- q_list[i]

  for (j in 1:5){
    TT <- T_list[j]

    for (k in 1:4){

      for (l in 2:2){

        cat("The current iteration is q:",q,"TT:",TT,"B_mat_type:",k,"deg_type:",l,"\n")
        load(paste0("sim_VHAR_q",q,"_T",TT,"_B_mat_type",k,"_deg_type",l,".RData"))

        iter_num <- length(VHAR_sim)

        # Kmeans:
        acc_kmeans <- round(1-mean(mapply(iter=1:iter_num,function(iter)
                              mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$kmeans[[m]]$err[1])))),3)
        ARI_kmeans <- round(mean(mapply(iter=1:iter_num,function(iter)
                            mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$kmeans[[m]]$err[2])))),3)

        # pam:
        acc_pam <- round(1-mean(mapply(iter=1:iter_num,function(iter)
                            mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$pam[[m]]$err[1])))),3)
        ARI_pam <- round(mean(mapply(iter=1:iter_num,function(iter)
                        mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$pam[[m]]$err[2])))),3)

        # noalpha:
        acc_noalpha <- round(1-mean(mapply(iter=1:iter_num,function(iter)
                                mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$noalpha[[m]]$err[1])))),3)
        ARI_noalpha <- round(mean(mapply(iter=1:iter_num,function(iter)
                              mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$noalpha[[m]]$err[2])))),3)
        
        # noalpha:
        acc_noalpha_pam <- round(1-mean(mapply(iter=1:iter_num,function(iter)
                              mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$noalpha_pam[[m]]$err[1])))),3)
        ARI_noalpha_pam <- round(mean(mapply(iter=1:iter_num,function(iter)
                              mean(mapply(m=1:4,function(m)VHAR_sim[[iter]]$noalpha_pam[[m]]$err[2])))),3)
        
        VHAR_summary <- rbind(VHAR_summary,data.frame(q=q,
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

write_xlsx(VHAR_summary, "VHAR_summary.xlsx")
