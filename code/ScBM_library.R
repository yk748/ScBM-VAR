#####################################################
# Phi generator
#####################################################
generate_directed_sbm <- function(q, sender, receiver, pref_matrix, selfTF=FALSE) {
  r <- nrow(pref_matrix)
  c <- ncol(pref_matrix)

  #------------- Check: n must be divisible by both r and c -------------#
  if (q %% r != 0 || q %% c != 0) {
    stop("n must be divisible by both nrow and ncol of the preference matrix.")
  }

  #------------- Heterogenous degrees -------------#
  y_deg <- rlnorm(q, meanlog = 2, sdlog = 1)
  z_deg <- rlnorm(q, meanlog = 2, sdlog = 1)
  
  #------------- For identifiability: -------------#
  for (k in 1:r){
    y_deg[which(sender==k)] <- y_deg[which(sender==k)]/sum(y_deg[which(sender==k)])
  }
  for (k in 1:c){
    z_deg[which(receiver==k)] <- z_deg[which(receiver==k)]/sum(z_deg[which(receiver==k)])
  }

  adj <- matrix(0, q, q)

  #------------- Sampling: -------------#
  p_mat <- matrix(NA,q,q)
  for (i in 1:q) {
    for (j in 1:q) {
      g_i <- sender[i]
      g_j <- receiver[j]
      p_mat[i,j] <- pref_matrix[g_i, g_j]*y_deg[i]*z_deg[j]
    }
  }
  p_mat <- p_mat/max(p_mat)
  adj <-1*(p_mat>matrix(runif(q^2,0,1),q,q))

  if(!selfTF){
    diag(adj) <- 0;
  }else{
    diag(adj) <- 1
  }

  return(list(adj = adj, sender_deg=y_deg, receiver_deg = z_deg))
}


Phi_generator <- function(q, B_list, s, phi, selfTF=FALSE,
                          n_comm, threshold=TRUE){

  Phi <- list()
  Adj <- list()
  w_Adj <- list()
  y_vec <- z_vec <- matrix(NA,q,s)
  for (m in 1:s){
    y_vec[,m] <- rep(seq(1,n_comm[2*(m-1)+1],by=1),each=q/n_comm[2*(m-1)+1])
    z_vec[,m] <- rep(seq(1,n_comm[2*m],by=1),each=q/n_comm[2*m])

    Adj[[m]] <- generate_directed_sbm(q, sender = y_vec[,m],
                                      receiver = z_vec[,m],
                                      pref_matrix=B_list[[m]], selfTF=selfTF)$adj;
    w <- matrix(runif(q^2,0.3,1),q,q)
    diag(w) <- 1
    w_Adj[[m]] <- w * Adj[[m]]

    #------------- Corrected not to divide by zero. It produced non-stationary errors -------------#
    In_deg <- rowSums(as.matrix(Adj[[m]])); # rowSums, O
    Out_deg <- colSums(as.matrix(Adj[[m]])); # colSums, P

    if (threshold){
      tau <- mean(Out_deg)
      In_deg <- In_deg + tau
      Out_deg <- Out_deg + tau
    }

    #------------- Set scaling factors: zero if degree is zero, else sqrt(1 / deg) -------------#
    out_scale <- ifelse(Out_deg == 0, 0, sqrt(1 / Out_deg))
    in_scale  <- ifelse(In_deg  == 0, 0, sqrt(1 / In_deg))

    #------------- Build scaling diagonal matrices -------------#
    D_out_inv_sqrt <- diag(out_scale, q);
    D_in_inv_sqrt  <- diag(in_scale, q)

    #------------- Final multiplication -------------#
    Phi_raw <- D_out_inv_sqrt %*% t(w_Adj[[m]]) %*% D_in_inv_sqrt
    Phi[[m]] <- phi[[m]]*(Phi_raw)
  }

  # Outputs
  output <- list();
  output$y_vec <- y_vec
  output$z_vec <- z_vec
  output$Adj <- Adj;
  output$W_Adj <- w_Adj
  output$Phi <- Phi;
  return(output)
}

#####################################################
# Power method for spectral radius
#####################################################
spectral_radius_power <- function(A, tol = 1e-6, max_iter = 1000) {
  x <- rnorm(ncol(A))
  x <- x / sqrt(sum(x^2))
  for (i in 1:max_iter) {
    x_new <- A %*% x
    lambda <- sqrt(sum(x_new^2))
    x_new <- x_new / lambda
    if (sqrt(sum((x_new - x)^2)) < tol) break
    x <- x_new
  }
  return(lambda)
}


#####################################################
# evaluate_clustering
#####################################################
evaluate_clustering <- function(true_labels, estimated_labels) {
  
  
  true_labels <- as.integer(factor(true_labels))
  estimated_labels <- as.integer(factor(estimated_labels))
  
  classes <- sort(unique(true_labels))
  k <- length(classes)
  
  perms <- permutations(n = k, r = k)
  
  best_error <- Inf
  best_ari <- NA
  best_labels <- NULL
  
  for (i in 1:nrow(perms)) {
    perm <- perms[i, ]
    
    #------------- Map estimated labels using current permutation -------------#
    label_map <- setNames(perm, classes)
    remapped <- sapply(estimated_labels, function(x) label_map[as.character(x)])
    
    #------------- Misclassification error -------------#
    err <- mean(remapped != true_labels)
    
    if (err < best_error) {
      best_error <- err
      best_labels <- remapped
      best_ari <- adjustedRandIndex(true_labels, remapped)
    }
  }
  
  return(list(
    err = c(best_error, best_ari),
    matched_labels = best_labels))
}


#####################################################
# Cross-validation for alpha
#####################################################
alpha_cv <- function(Phi_series,q,p,n_comm,type,clst,fold=5){
  
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
    
    # # ------------------------------- #
    # for (m in 1:s){
    #   tmp_L[,((m-1)*q+1):(m*q)] <- solve(sqrt(diag(rowSums(tmp_L[,((m-1)*q+1):(m*q)]),q) + diag(mean(colSums(tmp_L[,((m-1)*q+1):(m*q)])),q))) %*% tmp_L[,((m-1)*q+1):(m*q)] %*% solve(sqrt(diag(colSums(tmp_L[,((m-1)*q+1):(m*q)]),q) + diag(mean(colSums(tmp_L[,((m-1)*q+1):(m*q)])),q)))
    # }
    # # ------------------------------- #
    
    for(m in 1:s){
      Ky <- n_comm[(2*m-1)]; Kz <- n_comm[(2*m)]; Kyz <- min(Ky,Kz)
      svd_tmp <- svd((tmp_L[,((m-1)*q+1):(m*q)]),nu=Kyz,nv=Kyz)
      L_hat[,((m-1)*q+1):(m*q),l] <- svd_tmp$u[,1:Kyz] %*% diag(svd_tmp$d,Kyz) %*% t(Conj(svd_tmp$v[,1:Kyz]))
    }
  }
  
  #------------- Evaluate likelihood through PisCES -------------#
  # range of the candidate alphas
  #########################################################################################
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
              group_L[[1]] <- group_R[[s]] <- Blockbuster(V_aug,q,Kzy,clst)
              
            }else{
              group_L[[1]] <- Blockbuster(PisCES_la$VL_bar[[1]],q,Ky,clst)
              group_R[[s]] <- Blockbuster(PisCES_la$VR_bar[[s]],q,Kz,clst) 
            }
            
          }else{
            Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
            if (Ky == Kz){
              Kzy <- Ky
              V_aug <- rbind(PisCES_la$VR_bar[[(m-1)]],PisCES_la$VL_bar[[m]])
              group_L[[m]] <- group_R[[(m-1)]] <- Blockbuster(V_aug,q,Kzy,clst)
              
            }else{
              group_L[[m]] <- Blockbuster(PisCES_la$VL_bar[[m]],q,Ky,clst)
              group_R[[(m-1)]] <- Blockbuster(PisCES_la$VR_bar[[(m-1)]],q,Kz,clst)
              
            }
          }
        }
        
      }else if (type == "VHAR"){
        for (m in 1:s){
          if (m == 1){
            Ky <- n_comm[1]; Kz <- n_comm[(2*s)]
            group_L[[1]] <- Blockbuster(PisCES_la$VL_bar[[1]],q,Ky,clst)
            group_R[[s]] <- Blockbuster(PisCES_la$VR_bar[[s]],q,Kz,clst)
            
          }else{
            Ky <- n_comm[(2*m)-1]; Kz <- n_comm[2*(m-1)]
            if (Ky == Kz){
              Kzy <- Ky
              V_aug <- rbind(PisCES_la$VR_bar[[(m-1)]],PisCES_la$VL_bar[[m]])
              group_L[[m]] <- group_R[[(m-1)]] <- Blockbuster(V_aug,q,Kzy,clst)
              
            }else{
              group_L[[m]] <- Blockbuster(PisCES_la$VL_bar[[m]],q,Ky,clst)
              group_R[[(m-1)]] <- Blockbuster(PisCES_la$VR_bar[[(m-1)]],q,Kz,clst)
              
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
  M_large <- 1e8
  if (any(is.nan(loglik))){
    loglik[is.nan(loglik)] <- -M_large
  }
  
  loglik_final <- colSums(loglik)/fold
  alpha_hat <- alpha[which.max(loglik_final)]    
  
  output <- list()
  output$alpha_grid <- alpha
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
Blockbuster <- function(mat,q,K,clst=c("kmeans","pam","dbscan")){
  
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
    
    if (clst == "kmeans"){
      clut_result <- kmeans(X,centers=K)
      mat_1q <- cbind(1:q,clut_result$cluster)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }else if (clst == "pam"){
      clut_result <- pam(X,k=K)
      mat_1q <- cbind(1:q,clut_result$clustering)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }else { # clst == "dbscan"
      clut_result <- hdbscan(X,minPts = 2)
      mat_1q <- cbind(1:q,clut_result$cluster+1)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }
    
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
    
    if (clst == "kmeans"){
      clut_result <- kmeans(XX,centers=K)
      mat_1q <- cbind(1:q,clut_result$cluster)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }else if (clst == "pam"){
      clut_result <- pam(XX,k=K)
      mat_1q <- cbind(1:q,clut_result$clustering)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }else { # clst == "dbscan"
      clut_result <- hdbscan(XX,minPts=2)
      mat_1q <- cbind(1:q,clut_result$cluster+1)
      mat_re <- mat_1q[order(mat_1q[,2]),]
      
    }
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
PVAR_ols <- function(Yt,s,p=1){
  
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
