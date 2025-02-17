#####################################################
# Directed networks generator
#####################################################
DC_ScBM_generator <- function(q,s,n_comm,B_list,type,p,phi,tau=FALSE){
  
  #------------- Ym,Zm matrices -------------#
  y_vec <- z_vec <- array(NA,dim=c(q,s))
  if (type == "PVAR"){
    for (m in 1:s){
      if (m == 1){
        Ky <- n_comm[1]; Kz <- n_comm[(2*s)]
        if (Ky == Kz){
          Kyz <- Ky
          y_vec[,1] <- z_vec[,s] <- sample(rep(seq(1,Kyz,by=1),each=q/Kyz),
                                           size=q,replace=FALSE)
        }else{
          y_vec[,1] <- sample(rep(seq(1,Ky,by=1),each=q/Ky),
                              size=q,replace=FALSE)
          z_vec[,s] <- sample(rep(seq(1,Kz,by=1),each=q/Kz),
                              size=q,replace=FALSE)
        }

      }else{
        Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
        if (Ky == Kz){
          Kyz <- Ky
          y_vec[,m] <- z_vec[,(m-1)] <- sample(rep(seq(1,Kyz,by=1),each=q/Kyz),
                                               size=q,replace=FALSE)

        }else{
          y_vec[,m] <- sample(rep(seq(1,Ky,by=1),each=q/Ky),
                              size=q,replace=FALSE)
          z_vec[,(m-1)] <- sample(rep(seq(1,Kz,by=1),each=q/Kz),
                                  size=q,replace=FALSE)
        }
      }
    }
    
  }else if (type == "VHAR"){
    for (m in 1:s){
      if (m == 1){
        Ky <- n_comm[1]; Kz <- n_comm[2*s]
        y_vec[,1] <- sample(rep(seq(1,Ky,by=1),each=q/Ky),
                            size=q,replace=FALSE) 
        z_vec[,s] <- sample(rep(seq(1,Kz,by=1),each=q/Kz),
                            size=q,replace=FALSE)
        
      }else{
        Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
        if (Ky == Kz){
          Kyz <- Ky
          y_vec[,m] <- z_vec[,(m-1)] <- sample(rep(seq(1,Kyz,by=1),each=q/Kyz),
                                               size=q,replace=FALSE)
          
        }else{
          y_vec[,m] <- sample(rep(seq(1,Ky,by=1),each=q/Ky),
                              size=q,replace=FALSE)
          z_vec[,(m-1)] <- sample(rep(seq(1,Kz,by=1),each=q/Kz),
                                  size=q,replace=FALSE)
        }
      }
    }
  }
  
  # for (m in 1:s){
  #   y_vec[,m] <- sort(y_vec[,m])
  #   z_vec[,m] <- sort(z_vec[,m])
  # }
   
  Y <- Z <- list()
  for (m in 1:s){
    Ky <- n_comm[(2*m)-1]
    Y_mat <- array(0,dim=c(q,Ky))
    for (k in 1:Ky){
      Y_mat[y_vec[,m]==k,k] <- 1
    }

    Kz <- n_comm[2*m]
    Z_mat <- array(0,dim=c(q,Kz))
    for (k in 1:Kz){
      Z_mat[z_vec[,m]==k,k] <- 1
    }
    Y[[m]] <- Y_mat
    Z[[m]] <- Z_mat
  }
  
  #------------- Theta matrices -------------#
  Theta_O <- Theta_P <- list()
  for (m in 1:s){
    # From fastRG:
    y_deg <- rlnorm(q, meanlog = 2, sdlog = 1)
    z_deg <- rlnorm(q, meanlog = 2, sdlog = 1)
    
    # alpha <- 2.5
    # y_deg <- sample(alpha*(q/(1:q))^(1+alpha),q,replace=FALSE)
    # z_deg <- sample(rev(alpha*(q/(1:q))^(1+alpha)),q,replace=FALSE)
    
    # # Okay:
    # y_deg <- rep(1,q)
    # z_deg <- rev(rep(1,q))
    
    # y_deg <- seq(1,q)
    # z_deg <- rev(seq(1,q))  
    
    Ky <- n_comm[(2*m)-1]
    for (k in 1:Ky){
      y_deg[y_vec[,m]==k] <- y_deg[y_vec[,m]==k]/sum(y_deg[y_vec[,m]==k])
    }
    Theta_O[[m]] <- diag(y_deg)
    
    Kz <- n_comm[(2*m)]
    for (k in 1:Kz){
      z_deg[z_vec[,m]==k] <- z_deg[z_vec[,m]==k]/sum(z_deg[z_vec[,m]==k])
    }
    Theta_P[[m]] <- diag(z_deg)
  }
  
  #------------- Poisson sampling -------------#
  A <- list(); W <- list()
  for (m in 1:s){
    Ky <- n_comm[(2*m)-1]
    Kz <- n_comm[(2*m)]
    
    # B_mat <- diag(colSums(Theta_O[[m]] %*% Y[[m]])) %*% B_list[[m]] %*% diag(colSums(Theta_P[[m]] %*% Z[[m]]))
    B_mat <- t(Theta_O[[m]] %*% Y[[m]]) %*% (Theta_O[[m]] %*% Y[[m]]) %*% B_list[[m]] %*% t(Theta_P[[m]] %*% Z[[m]]) %*% (Theta_P[[m]] %*% Z[[m]])
    Prob_B_mat <- B_mat/sum(B_mat)
    ref_B <- matrix(c(1:length(B_mat)),Ky,Kz)
    
    Y_sum <- colSums(Theta_O[[m]] %*% Y[[m]]); Z_sum <- colSums(Theta_P[[m]] %*% Z[[m]])
    Y_tilde <- array(0,dim(Y[[m]])); Z_tilde <- array(0,dim(Z[[m]]))
    for (j in 1:Ky){
      if (Y_sum[j]!=0){
        Y_tilde[,j] <- (Theta_O[[m]] %*% Y[[m]])[,j]/Y_sum[j]
      }
    }
    for (j in 1:Kz){
      if (Z_sum[j]!=0){
        Z_tilde[,j] <- (Theta_P[[m]] %*% Z[[m]])[,j]/Z_sum[j]
      }
    }
    
    dnst <- 1
    num_ent <- dnst*q^2
    A[[m]] <- array(0,dim=c(q,q))
    add_ref <- sample(c(1:(Ky*Kz)),num_ent,replace=TRUE,prob=as.vector(Prob_B_mat))
    for (i in 1:num_ent){
      coord_B <- which(ref_B==add_ref[i],arr.ind=TRUE)
      coord_y <- sample(1:q,1,replace=TRUE,prob=Y_tilde[,coord_B[1]])
      coord_z <- sample(1:q,1,replace=TRUE,prob=Z_tilde[,coord_B[2]])
      
      A[[m]][coord_y,coord_z] <- A[[m]][coord_y,coord_z] + 1
    }
    A[[m]] <- 1*(A[[m]]!=0)
    # diag(A[[m]]) <- rep(1,q)
    W[[m]] <- A[[m]]*runif(q^2,0.3,1)
    # diag(W[[m]]) <- rep(1,q)
  }
  
  #------------- Laplacian matrices & VAR Transition matrices -------------#
  L <- list()
  Phi <- list()
  for (m in 1:s){
    
    if (isFALSE(tau)){
      tau <- 0
    }else{
      tau <- mean(colSums(W[[m]]))
    }
    
    O_tau <- colSums(W[[m]]) + rep(tau,q)
    P_tau <- rowSums(W[[m]]) + rep(tau,q)
    O_sqrt_inv <- matrix(0,nrow=q,ncol=q)
    P_sqrt_inv <- matrix(0,nrow=q,ncol=q)
    for (i in 1:q){
      if(P_tau[i] != 0){
        P_sqrt_inv[i,i] <- 1/sqrt(P_tau[i])
      }else{
        P_sqrt_inv[i,i] <- 0 
      }
      
      if(O_tau[i] != 0){
        O_sqrt_inv[i,i] <- 1/sqrt(O_tau[i])
      }else{
        O_sqrt_inv[i,i] <- 0 
      }
    }
    
    Phi_m <- array(NA,dim=c(q,(q*p)))
    L_m <- matrix(0,q,q)
    for (h in 1:p){
      Phi_m[,((h-1)*q+1):(h*q)] <- phi[[m]][h]* t(O_sqrt_inv %*% W[[m]] %*% P_sqrt_inv)
      L_m <- L_m + t(Phi_m[,((h-1)*q+1):(h*q)])/phi[[m]][h]
    }
    Phi[[m]] <- Phi_m
    L[[m]] <- L_m
  }
  
  #------------- Outputs -------------#
  output <- list()
  output$y_vec <- y_vec
  output$z_vec <- z_vec
  output$Y <- Y
  output$Z <- Z
  output$Theta_O <- Theta_O
  output$Theta_P <- Theta_P
  output$W <- W
  output$A <- A
  output$L <- L
  output$Phi <- Phi
  return(output)
}


#####################################################
# Cross-validation for alpha
#####################################################
alpha_cv <- function(Phi_series,q,p,n_comm=n_comm,type,fold=5){
  
  #------------- Construct completed adjacency matrices -------------#
  # Divide q (diagonal) + comb(q,2) pairs into S fold
  Id_diag <- mapply(i=1:q,FUN=function(i)q*(i-1)+i)
  qc2 <- q*(q-1)/2
  
  s <- length(n_comm)/2 # # of seasons
  Idx <- array(NA,dim=c(s,qc2,fold))
  for (l in 1:fold){
    for (m in 1:s){
      Idx[m,,l] <- sample(setdiff(c(1:q^2),Id_diag),qc2, replace = FALSE, prob = NULL)
    }
  }
  
  # Complete Laplacian matrices by SVD
  L_hat <- array(NA,dim=c(q,(s*q),fold));
  for(l in 1:fold){
    tmp_L <- matrix(mapply(x=1:s,function(x) 
      selectPairMatrix(q,Idx[x,,l],Phi_series[,((x-1)*q+1):(x*q)])),
      nrow=q,ncol=(s*q))
    
    for(m in 1:s){
      Ky <- n_comm[(2*m-1)]; Kz <- n_comm[(2*m)]; Kyz <- min(Ky,Kz)
      svd_tmp <- svd((tmp_L[,((m-1)*q+1):(m*q)]),nu=Kyz,nv=Kyz)
      L_hat[,((m-1)*q+1):(m*q),l] <- svd_tmp$u[,1:Kyz] %*% diag(svd_tmp$d,Kyz) %*% t(Conj(svd_tmp$v[,1:Kyz]))/0.9
    }
  }
  
  #------------- Evaluate likelihood through PisCES -------------#
  # range of the candidate alphas
  #########################################################################################
  # This has to be investigated:
  alpha_max <- 1/(4*sqrt(2)+2)
  alpha <- exp(seq(log(0.01*alpha_max),
                   log(alpha_max),length.out=20))
  #########################################################################################
  loglik <- array(NA,dim=c(fold,length(alpha)))
  
  for(l in 1:fold){
    L_hat_l <- L_hat[,,l]
    
    for (a in 1:length(alpha)){
      # For a fixed alpha and fold, run PisCES
      PisCES_la <- PisCES(L_hat_l,q,p,n_comm,s,alpha[a])
      
      #------------- reconstruct Laplacian -------------#
      loglik_la <- vector("numeric",s)
      
      # Identify co-group structures
      group_L <- list(); group_R <- list();
      if (type == "PVAR"){
        for (m in 1:s){
          if (m == 1){
            Ky <- n_comm[1]; Kz <- n_comm[(2*s)]
            if (Ky == Kz){
              Kzy <- Ky
              V_aug <- rbind(PisCES_la$VR_bar[[s]],PisCES_la$VL_bar[[1]])
              group_L[[1]] <- group_R[[s]] <- Blockbuster(V_aug,q,Kzy)
              
            }else{
              group_L[[1]] <- Blockbuster(PisCES_la$VL_bar[[1]],q,Ky)
              group_R[[s]] <- Blockbuster(PisCES_la$VR_bar[[s]],q,Kz) 
            }

          }else{
            Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
            if (Ky == Kz){
              Kzy <- Ky
              V_aug <- rbind(PisCES_la$VR_bar[[(m-1)]],PisCES_la$VL_bar[[m]])
              group_L[[m]] <- group_R[[(m-1)]] <- Blockbuster(V_aug,q,Kzy)

            }else{
              group_L[[m]] <- Blockbuster(PisCES_la$VL_bar[[m]],q,Ky)
              group_R[[(m-1)]] <- Blockbuster(PisCES_la$VR_bar[[(m-1)]],q,Kz)

            }
          }
        }
        
      }else if (type == "VHAR"){
        for (m in 1:s){
          if (m == 1){
            Ky <- n_comm[1]; Kz <- n_comm[(2*s)]
            group_L[[1]] <- Blockbuster(PisCES_la$VL_bar[[1]],q,Ky)
            group_R[[s]] <- Blockbuster(PisCES_la$VR_bar[[s]],q,Kz)
            
          }else{
            Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
            if (Ky == Kz){
              Kzy <- Ky
              V_aug <- rbind(PisCES_la$VR_bar[[(m-1)]],PisCES_la$VL_bar[[m]])
              group_L[[m]] <- group_R[[(m-1)]] <- Blockbuster(V_aug,q,Kzy)
              
            }else{
              group_L[[m]] <- Blockbuster(PisCES_la$VL_bar[[m]],q,Ky)
              group_R[[(m-1)]] <- Blockbuster(PisCES_la$VR_bar[[(m-1)]],q,Kz)
              
            }
          }
        }
        
      }
      
      # Reconstruct estimated Laplacian matrix
      L_tilde_la <- array(NA,dim=c(q,q,s))
      for (m in 1:s){
        # Compute node-community incidence matrix
        num_Y <- length(unique(group_L[[m]]$group_1q))
        Y_mat <- matrix(0,q,num_Y)
        for (g in 1:num_Y){
          Y_mat[,g] <- 1*(group_L[[m]]$group_1q == g)
        }
        num_Z <- length(unique(group_R[[m]]$group_1q))
        Z_mat <- matrix(0,q,num_Z)
        for (g in 1:num_Z){
          Z_mat[,g] <- 1*(group_R[[m]]$group_1q == g)
        }
        
        # Compute in- and out-degrees & re-order them
        L_hat_lam <- L_hat_l[,((m-1)*q+1):(m*q)]
        d_y <- rowSums(L_hat_lam)
        d_z <- colSums(L_hat_lam)
        
        # Node-community probability
        B_yz <- matrix(NA,num_Y,num_Z)
        for (y in 1:num_Y){
          for (z in 1:num_Z){ 
            n_yz <- sum(d_y[group_L[[m]]$group_1q==y] %*% t(d_z[group_R[[m]]$group_1q==z]))
            B_yz[y,z] <- sum(L_hat_lam[group_L[[m]]$group_1q==y,
                                       group_R[[m]]$group_1q==z])/n_yz
          }
        }
        
        # Compute probability matrix
        L_tilde_la[,,m] <- diag(d_y,q) %*% Y_mat %*% B_yz %*% t(Z_mat) %*% diag(d_z,q)
        
        # Compute 2nd-order approximation of von Neumann entropy
        loglik_la[m] <- sum(diag(Phi_series[,((m-1)*q+1):(m*q)]))/q*(1 -sum(diag(L_tilde_la[,,m]))/q)
      }
      loglik[l,a] <- sum(loglik_la)
    }
  }
  
  #------------- Choose alpha -------------#
  loglik_final <- colSums(loglik)/fold
  alpha_hat <- alpha[which.max(loglik_final)]    
  
  output <- list()
  output$alpha_hat <- alpha_hat
  output$likelihood <- loglik_final
  return(output)
}

#####################################################
# Matrix selector
#####################################################
selectPairMatrix = function(q, index, mat){
  
  y <- rep(1,q^2)
  y[index] <- 0
  A0 <- matrix(y,nrow=q,ncol=q)
  # Note: A should be from (row) -> to (column)
  # But mat is a transition matrix (to (row) <- from (column)) 
  # So it should be transposed.
  A <- t(mat) * A0
  
  return(A)
}

#####################################################
# Blockbuster: Spectral clustering
#####################################################
Blockbuster <- function(mat,q,K){
  
  qq <- dim(mat)[1]
  
  # Stack column if needed. 
  # Otherwise, proceed separately
  if (qq == q){
    X <- mat
    sd_X <- mapply(x=1:q,function(x)norm(X[x,],"2"))
    for (i in 1:q){
      if(sd_X[i]!=0){
        X[i,] <- X[i,]/sd_X[i]
      } 
    }
    
    k_means <- kmeans(X,centers=K)
    mat_1q <- cbind(1:q,k_means$cluster)
    mat_re <- mat_1q[order(mat_1q[,2]),]
    
  }else{
    X1 <- mat[1:q,]; X2 <- mat[(q+1):(2*q),]
    sd_X1 <- mapply(x=1:q,function(x)norm(X1[x,],"2"))
    sd_X2 <- mapply(x=1:q,function(x)norm(X2[x,],"2"))
    for (i in 1:q){
      if(sd_X1[i]!=0){
        X1[i,] <- X1[i,]/sd_X1[i]
      } 
      if(sd_X2[i]!=0){
        X2[i,] <- X2[i,]/sd_X2[i]
      } 
    }
    XX <- cbind(X1,X2)
    
    k_means <- kmeans(XX,centers=K)
    mat_1q <- cbind(1:q,k_means$cluster)
    mat_re <- mat_1q[order(mat_1q[,2]),]
  }
  
  output <- list()
  output$subj <- mat_re[,1]
  output$group_re <- mat_re[,2]
  output$group_1q <- mat_1q[,2]
  return(output)
}

#####################################################
# PisCES: Singular vectors smoother
#####################################################
PisCES <- function(Phi,q,p,n_comm,s,alpha){
  
  #------------- Initialization -------------#
  Laplace <- array(0,dim=c(q,q,s));
  UL_hat <- array(NA,dim=c(q,q,s)); UR_hat <- array(NA,dim=c(q,q,s))
  SVD_K <- list()
  
  for (m in 1:s){
    for (l in 1:p){
      Laplace[,,m] <- Laplace[,,m] + t(Phi[((l-1)*q+1):(l*q),((m-1)*q+1):(m*q)]) 
    }
    Ky <- n_comm[(2*m-1)]; Kz <- n_comm[(2*m)]
    SVD_K[[m]] <- svd(Laplace[,,m])
    UL_hat[,,m] <- SVD_K[[m]]$u[,1:Ky] %*% t(SVD_K[[m]]$u[,1:Ky])
    UR_hat[,,m] <- SVD_K[[m]]$v[,1:Kz] %*% t(SVD_K[[m]]$v[,1:Kz])
  }
  
  #------------- Left singular vector smoothing -------------#
  UL_bar_cur <- UL_hat; UL_bar_next <- array(0,dim=dim(UL_hat)); svd_L <- list()
  eps <- 1e-5; max_iter <- 1000; iter_L <- 0
  while (iter_L <= max_iter){
    iter_L <- iter_L + 1
    
    diff <- 0
    for (m in 1:s){
      Ky <- n_comm[(2*m-1)]
      
      if (m == 1){
        svd_L[[1]] <- svd(UL_hat[,,1] + alpha*UL_bar_cur[,,2],nu=Ky,nv=Ky)
        UL_bar_next[,,1] <- svd_L[[1]]$u %*% t(svd_L[[1]]$u)
        
      }else if (m == s){
        svd_L[[s]] <- svd(alpha*UL_bar_cur[,,(s-1)] + UL_hat[,,s],nu=Ky,nv=Ky)
        UL_bar_next[,,s] <- svd_L[[s]]$u %*% t(svd_L[[s]]$u)
        
      }else{
        svd_L[[m]] <- svd(alpha*UL_bar_cur[,,(m-1)] 
                          + UL_hat[,,m] +alpha*UL_bar_cur[,,(m+1)],nu=Ky,nv=Ky)
        UL_bar_next[,,m] <- svd_L[[m]]$u %*% t(svd_L[[m]]$u)
        
      }
      diff <- diff + norm(UL_bar_next[,,m] - UL_bar_cur[,,m],"F")
    }
    
    if (diff < eps){
      break
    }else{
      UL_bar_cur <- UL_bar_next
    }
  }
  
  #------------- Right singular vector smoothing -------------#
  UR_bar_cur <- UR_hat; UR_bar_next <- array(0,dim=dim(UR_hat)); svd_R <- list()
  eps <- 1e-5; max_iter <- 1000; iter_R <- 0
  while (iter_R <= max_iter){
    iter_R <- iter_R + 1
    
    diff <- 0
    for (m in 1:s){
      Kz <- n_comm[(2*m)]
      
      if (m == 1){
        svd_R[[1]] <- svd(UR_hat[,,1] + alpha*UR_bar_cur[,,2],nu=Kz,nv=Kz)
        UR_bar_next[,,1] <- svd_R[[1]]$v %*% t(svd_R[[1]]$v)
        
      }else if (m==s){
        svd_R[[s]] <- svd(alpha*UR_bar_cur[,,(s-1)] + UR_hat[,,s],nu=Kz,nv=Kz)
        UR_bar_next[,,s] <- svd_R[[s]]$v %*% t(svd_R[[s]]$v)
        
      }else{
        svd_R[[m]] <- svd(alpha*UR_bar_cur[,,(m-1)] 
                          + UR_hat[,,m] +alpha*UR_bar_cur[,,(m+1)],nu=Kz,nv=Kz)
        UR_bar_next[,,m] <- svd_R[[m]]$v %*% t(svd_R[[m]]$v)
        
      }
      diff <- diff + norm(UR_bar_next[,,m] - UR_bar_cur[,,m],"F")
    }
    
    if (diff < eps){
      break
    }else{
      UR_bar_cur <- UR_bar_next
    }
  }
  
  VL_bar <- list(); VR_bar <- list()
  for (m in 1:s){
    Ky <- n_comm[(2*m-1)]; Kz <- n_comm[(2*m)]
    
    svd_UL <- svd(UL_bar_next[,,m],nu=Ky,nv=Ky)
    svd_UR <- svd(UR_bar_next[,,m],nu=Kz,nv=Kz)
    
    VL_bar[[m]] <- svd_UL$u[,1:Ky] %*% sqrt(diag(svd_UL$d[1:Ky],Ky))
    VR_bar[[m]] <- svd_UR$v[,1:Kz] %*% sqrt(diag(svd_UR$d[1:Kz],Kz))
  }
  
  #------------- Report output -------------#
  output <- list()
  output$UL <- UL_bar_next
  output$UR <- UR_bar_next
  output$iter_L <- iter_L
  output$iter_R <- iter_R
  output$VL_bar <- VL_bar
  output$VR_bar <- VR_bar
  return(output)
}

##############################################################
# VAR(p)_generator
##############################################################
VAR_sim <- function(TT, A, p, Sigma){
  
  q <- dim(A)[1];
  burn <- 500;
  
  inno <- mvrnorm(n=TT+burn, rep(0, q), Sigma);
  init <- mvrnorm(n=p, rep(0, q), Sigma);
  init <- matrix(init, nrow=p);
  
  #------------- Find index for previous observations -------------#
  j <- 1;
  # ar term
  id <- seq(from= j+p-1, to = j, by=-1);
  
  Y <- matrix(0, (TT+burn), q);
  for(t in 1:(TT+burn)){
    Y[t,] = A%*%as.vector(t(init[id,])) + inno[t,];
    init = rbind(init[-1,], Y[t,]);
  }
  
  return(t(Y[-(1:burn),])) # Final data is k*T matrix
}


#####################################################
# PVAR_generator
#####################################################
PVAR_generator = function(TT, s, p, Phi, Sigma){
  
  q <- dim(Phi)[1];
  burn <- 500;
  
  innovations <- mvrnorm(TT+burn, rep(0,q), Sigma);
  
  Yt <- matrix(0,TT+burn,q);
  Yt[1:p,] <- innovations[1:p,];
  
  # for (t in 2:(TT+burn)){
  #   m <- (t-1) %% s + 1
  #   Yt[t,] <- Phi[,((m-1)*q+1):(m*q)] %*% Yt[(t-1),] + innovations[t,]
  # }

  Theta <- phitoTheta(Phi,s,q,p)
  for(t in (p+1):(TT+burn)){

    j <- t %% s;
    j <- ifelse(j == 0, s, j);
    idx_y <- seq(from= t-1, to = t-p, by=-1);
    vec_yt <- as.vector(Yt[idx_y,]);

    for(k in 1:q){
      idx_q = seq(from=(k-1)*(p*q)+1, to=k*(p*q), by=1);
      # cat(idx_q,"\n")
      Yt[t,k] = Theta[j,idx_q] %*% vec_yt + innovations[t,k];
    }
  }
  
  Yt <- t(Yt[-(1:burn),]);
  return(Yt)
}

#####################################################
# phitoTheta
#####################################################
phitoTheta = function(Phi,s,q,p){
  
  Theta <- matrix(0, s, q^2*p);
  
  for(i in 1:s){
    idx <- seq(from=(i-1)*q+1, to= i*q, by=1);
    Theta[i,] <- as.vector(Phi[,idx]);
  }
  
  return(Theta)
}


#####################################################
# VHAR_generator
#####################################################
VHAR_generator = function(TT, Phi, Sigma){
  
  q <- dim(Phi)[1];
  burn <- 500;
  
  Daily <- Phi[,(1:q)];
  Weekly <- Phi[,((q+1):(2*q))];
  Monthly <- Phi[,((2*q+1):(3*q))];
  
  Phi_aug = Daily + 1/5*Weekly + 1/22*Monthly;
  for(i in 1:4){
    Phi_aug <- cbind(Phi_aug, 1/5*Weekly + 1/22*Monthly)
  }
  for(i in 1:17){
    Phi_aug <- cbind(Phi_aug, 1/22*Monthly)
  }
  
  innovations <- mvrnorm(TT+burn, rep(0,q), Sigma);
  initial <- mvrnorm(22, rep(0,q), Sigma);
  
  Phi_idx <- seq(from= 22, to = 1, by=-1);
  
  Yt <- matrix(0, (TT+burn),q);
  for(t in 1:(TT+burn)){
    Yt[t,] <- Phi_aug %*% as.vector(t(initial[Phi_idx,])) + innovations[t,];
    initial <- rbind(initial[-1,], Yt[t,]);
  }
  
  Yt <- t(Yt[-(1:burn),]);
  return(Yt) # Final data is dim*T matrix
}


#####################################################
# PVAR estimation through OLS:
#####################################################
PVAR_ols_new <- function(Yt,s,p=1){
  
  TT <- ncol(Yt);
  q <- nrow(Yt);
  N <- TT/s
  
  Phi_hat <- array(NA,dim=c(q,q*s))
  for (m in 1:s){
    idx <- which((c(1:TT) - 1) %% s + 1 == m)[-1]
    Z <- Yt[,idx]
    X <- Yt[,(idx-1)]
    
    eig_val <- eigen(X%*%t(X))$values
    min_idx <- which(eig_val < 0)
    if(length(min_idx)==0){
      eps <- 0
    }else{
      # eps <- 50*abs(eig_val[max(min_idx)])
      eps <- 1e-8
    }
    
    Phi_hat[,(q*(m-1)+1):(q*m)] <- solve(X %*% t(X) + eps*diag(1,q)) %*% X %*% t(Z)
  }
  return(list(Phi_hat=Phi_hat))
}


PVAR_ols = function(Yt, s, p){
  
  TT <- ncol(Yt);
  q <- nrow(Yt);
  if(is.null(q)){
    q <- 1
  }
  
  # st <- (floor(p/s-10e-6)+1)*s; 
  st <- 1
  hat <- numeric(q^2*p);
  Phi_hat <- NULL;
  
  m = floor(TT/s);
  for(i in 1:s){
    
    id <- seq(from=st+i, to=TT, by=s);
    Y1 <- as.vector(Yt[,id]);
    len_id <- length(id);
    X1 <- matrix(0, p*q, len_id); 
    
    for(j in 1:len_id){
      
      id_ar <- seq(from=id[j]-1 , to = id[j]-p ,by=-1);
      # cat(id_ar,"\n")
      X1[,j] <- as.vector(Yt[,id_ar]);
    }
    
    # eps <- 0
    eig_val <- eigen(X1%*%t(X1))$values
    min_idx <- which(eig_val < 0)
    if(length(min_idx)==0){
      eps <- 0
    }else{
      # eps <- 50*abs(eig_val[max(min_idx)])
      eps <- 1e-8
    }
    
    # Replace the line using kronecker with:
    hat <- numeric(q*p*q)
    for (i in 1:q) {
      start_idx <- (i-1)*p*q + 1
      end_idx <- i*p*q
      # cat(start_idx,end_idx,"\n")
      hat[start_idx:end_idx] <- solve(X1 %*% t(X1) + eps * diag(1, q)) %*% X1 %*% Y1[(i-1)*len_id + 1:len_id]
    }
    
    a   <- matrix(hat, nrow=q)
    Phi_hat <- rbind(Phi_hat, a);
  }
  init <- phitoTheta(t(Phi_hat),s,q,p);
  
  return(list(Phi_hat=t(Phi_hat), init=init) );
}

#####################################################
# Periodic mean estimation in PVAR:
#####################################################
PCstat = function(y,s){
  
  # y is 1*n vector
  q <- floor(length(y)/s);
  y1 <- y[1:(q*s)];
  z <- matrix(y1, s);
  
  PCstd <- apply(z, 1, sd);
  PCmean <- apply(z, 1, mean);
  
  N1 <- dim(z)[2];
  up.mean <- qt(.975, df=N1-1)*PCstd/sqrt(N1);
  up.sd <- sqrt((N1-1)*PCstd^2/qchisq(.025, N1-1));
  dn.sd <- sqrt((N1-1)*PCstd^2/qchisq(.975, N1-1));
  
  return(list(mean=PCmean, sd = PCstd, 
              up.mean=PCmean+up.mean, dn.mean=PCmean-up.mean, 
              up.sd=up.sd, dn.sd=dn.sd));
}


#####################################################
# VHAR estimation through OLS:
#####################################################
VHAR_ols <- function(Yt){
  
  Yt_tr <- t(Yt);
  TT <- dim(Yt_tr)[1]
  q <- dim(Yt_tr)[2]
  
  X_cal <- VHAR_X_matrix(Yt_tr,TT,q)
  y_cal <- VHAR_y_matrix(Yt_tr)
  
  eig_val <- eigen(X_cal %*% t(X_cal))$values
  min_idx <- which(eig_val < 0)
  if(length(min_idx)==0){
    eps <- 0
  }else{
    # eps <- 50*abs(eig_val[max(min_idx)])
    eps <- 1e-8
  }
  
  Phi_hat <- y_cal %*% t(X_cal) %*% solve(X_cal %*% t(X_cal) + eps*diag(1,(3*q)));
  return(list(Phi_hat = Phi_hat, 
              D_hat = Phi_hat[,(1:q)], W_hat = Phi_hat[,((q+1):(2*q))], M_hat = Phi_hat[,-(1:(2*q))]))
}


#####################################################
# X matrix construction for VHAR estimation:
#####################################################
VHAR_X_matrix <- function(Yt_tr,TT,q){
  D <- t(Yt_tr)[,-(1:21)]
  D <- D[,-(ncol(D))]
  
  W <- matrix(0, nrow = TT-22, ncol=q)
  W.data <- Yt_tr[-(1:17),]
  for(i in 1:ncol(W)){
    temp <- na.omit( stats::filter(W.data[,i], rep(1/5, 5), sides=1) )
    W[,i] <- temp[-length(temp)]
  }
  W <- t(W)
  
  M <- matrix(0, nrow = TT-22, ncol=q)
  M.data <- Yt_tr
  for(i in 1:ncol(M)){
    temp <- na.omit(stats::filter(M.data[,i], rep(1/22,22), sides=1) )
    M[,i] <- temp[-length(temp)]
  }
  M <- t(M)
  
  return(rbind(D,W,M))
}

#####################################################
# Y matrix construction for VHAR estimation:
#####################################################
VHAR_y_matrix <- function(Yt_tr){
  return( t(Yt_tr)[,-(1:22)] )
}

#####################################################
# Sigma estimation in VHAR:
#####################################################
VHAR_sigma = function(y, x, hA1){
  
  Resi <- y - hA1%*%x;
  Sigma_z <- Resi%*%t(Resi)/ncol(Resi);
  out <- list();
  out$Resi <- Resi;
  out$Sigma_z <- Sigma_z;
  
  return(out);
}

#####################################################
# VAR estimation through OLS:
#####################################################
VAR_OLS = function(Yt, p){
  
  Yt <- Yt - rowMeans(Yt);
  q = dim(Yt)[1];
  TT = dim(Yt)[2];
  T1 = TT-p;
  
  X1 <- matrix(0,q*p, T1); 
  Y1 <- matrix(0,q, T1);
  for(j in 1:T1){
    
    id <- seq(from= j+p-1, to = j, by=-1);
    x <- as.vector(Yt[,id]);
    X1[,j] <- x;
    Y1[,j] <- Yt[,(j+p)];
  }
  
  Phi_hat <- Y1%*%t(X1)%*%solve(X1%*%t(X1));
  
  Resi <- matrix(0,q, T1); 
  for(j in 1:T1){
    
    id <- seq(from= j+p-1, to = j, by=-1);
    x <- as.vector(Yt[,id]);
    Resi[,j] =  Yt[,p+j]  - Phi_hat%*%x;
  }
  Sigma_z <- Resi %*% t(Resi)/T1;
  bic <- T1*log(det(Sigma_z)) + log(T1)* sum(Phi_hat != 0);
  
  return(list(Phi_hat = Phi_hat, Sigma=Sigma_z, bic=bic, p=p))
}


#####################################################
# Sigma estimation in VAR:
#####################################################
VAR_sigma = function(y, hA1){
  k = dim(y)[1];
  T = dim(y)[2];
  p = dim(hA1)[2]/k;
  T1 = T - p;
  
  # Calculate residuals
  Resi = matrix(0, k, T1);
  for(j in 1:T1){
    # ar term
    id = seq(from= j+p-1, to = j, by=-1);
    Resi[,j] =  y[,p+j]  - hA1%*%as.vector(y[,id]);
  }
  
  # Estimate Sigma
  Sigma_z = Resi%*%t(Resi)/T;
  out = list();
  out$Resi = Resi;
  out$Sigma_z = Sigma_z;
  return(out);
}


#####################################################
# OLS estimation:
#####################################################
const_ols = function(x,y, A.lasso){ 
  J <- 1*(A.lasso != 0); p=3;
  vecA1 <- c(J);
  A1 <- matrix(vecA1, nrow=q);
  k2 <- sum(vecA1);
  R <- matrix(0, (q^2)*p, k2);
  
  for(i in 1:(q^2*p)){
    if(vecA1[i] == 1){ 
      id <- sum(vecA1[1:i]);  
      R[i,id] <- 1;
    }  
  }
  
  Sigma_update <- VHAR.sigma(y, x, A.lasso)$Sigma_z;
  sv <- svd(Sigma_update)
  szInv <- sv$u%*%diag(1/sv$d)%*%t(sv$v)
  varA1 <- solve(t(R) %*% kronecker(x %*% t(x), szInv) %*% R);
  A1.db <- matrix( R %*% varA1 %*% t(R) %*% kronecker(x, szInv) %*% as.vector(y), nrow=q);
  return(A1.db)
}


#####################################################
# Calculate NMI:
#####################################################
calculate_nmi <- function(labels1, labels2) {
  contingency_table <- table(labels1, labels2)
  
  # Calculate Mutual Information (MI)
  joint_prob <- contingency_table / sum(contingency_table)
  row_marginals <- rowSums(joint_prob)
  col_marginals <- colSums(joint_prob)
  
  mi <- sum(
    joint_prob[joint_prob > 0] * 
      log(joint_prob[joint_prob > 0] / (row_marginals %*% t(col_marginals))[joint_prob > 0])
  )
  
  # Calculate entropy
  entropy <- function(prob) {
    -sum(prob[prob > 0] * log(prob[prob > 0]))
  }
  
  h1 <- entropy(row_marginals)
  h2 <- entropy(col_marginals)
  
  # Calculate NMI
  #  nmi <- mi / sqrt(h1 * h2)
  nmi <- 2*mi / (h1 + h2)
  return(nmi)
}


#####################################################
# Refiner
#####################################################
Refiner <- function(result_list){
  
  # Check the number of columns in the label matrix:
  ncol <- length(result_list)
  mat_refined <- array(NA,dim=c(q,ncol))
  
  # Initialization:
  mat_refined[,1] <- result_list[[1]]$group_1q
  
  # Find rearrangement:
  for (j in 2:ncol){
    num_comm <- length(unique(result_list[[j]]$group_1q))
    entry <- c(1:num_comm)
    perm_mat <- matrix(unlist(permn(entry)),ncol=num_comm,byrow=TRUE)
    
    submat_refined <- array(NA,dim=c(q,dim(perm_mat)[1]))
    for (i in 1:dim(perm_mat)[1]){
      mat_tmp <- vector("numeric",q)
      for (k in 1:length(entry)){
        mat_tmp[which(result_list[[j]]$group_1q == entry[k])] <- perm_mat[i,k]
      }
      submat_refined[,i] <- mat_tmp
    }
    
    # # Criteria: Number of equivalent labels 
    # equiv <- (submat_refined == mat_refined[,(j-1)])
    # sum_equiv <- colSums(equiv)
    # mat_refined[,j] <- submat_refined[,which.min(sum_equiv)]
    
    # Criteria: Minimizing absolute difference
    MAD <- abs(submat_refined - mat_refined[,(j-1)])
    sum_MAD <- colSums(MAD)
    mat_refined[,j] <- submat_refined[,which.min(sum_MAD)]
    
    result_list[[j]]$group_1q <- mat_refined[,j]
    mat_1q <- cbind(1:q,mat_refined[,j])
    mat_re <- mat_1q[order(mat_1q[,2]),]
    result_list[[j]]$group_re <- mat_re[,2]
    result_list[[j]]$subj <- mat_re[,1]
  }
  
  # Return output:
  return(result_list)
}

#####################################################
# Companion matrix from VAR(p) to make VAR(1)
#####################################################
companion_form_phi = function(Psi,d,p){
  
  comp_Psi <- array(0,dim=c((p*d),(p*d)))
  if (p == 1){
    comp_Psi = Psi[,,1]
  }
  else{
    for (i in 1:p){
      for (j in 1:p){
        if(i == 1){
          comp_Psi[1:d,((j-1)*d+1):(j*d)] <- Psi[,,j]
        }
        else if ( (i-1) == j & i != 1 ){
          comp_Psi[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- diag(1,d)
        }
      }
    }
  }
  
  return(comp_Psi)
}
