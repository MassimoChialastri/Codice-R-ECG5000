FKM.mahal <- function(X, k, m = 2, stand = 0, startU = NULL,
                      conv = 1e-9, maxit = 1e6, clus_var){
  library(stats)
  # scale data if stand = 1
  if(stand == 1){
    X <- scale(X)
  }
  # initialize fuzzy partition matrix
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(startU)){
    U <- matrix(runif(n * k), nrow = n, ncol = k)
    U <- U / rowSums(U)
  }else{
    U <- startU
  }
  # iterate until convergence or maximum iterations
  for(iter in 1:maxit){
    # update cluster centers
    H <- matrix(0, nrow = k, ncol = d)
    for(j in 1:k){
      for(i in 1:n){
        H[j, ] <- H[j, ] + (U[i, j]^m) * X[i, ]
      }
      H[j, ] <- H[j, ] / sum(U[, j]^m)
    }
    # update fuzzy partition matrix
    dist_matrix <- matrix(0, n, k)
    for(j in 1:k){
      dist_matrix[,j] <- sqrt(mahalanobis(X, center = H[j,], 
                                          cov = clus_var[[j]]))
    }
    U_new <- matrix(0, nrow = n, ncol = k)
    for(i in 1:n){
      for(j in 1:k){
        denom <- sum((dist_matrix[i, j] / dist_matrix[i, ])^(2/(m-1)))
        U_new[i, j] <- 1 / denom
      }
    }
    # check for convergence
    if(sum(abs(U_new - U)) < conv) {
      break
    }
    U <- U_new
  }
  out <- list()
  out$U <- U # membership degree matrix
  out$H <- H # prototype matrix
  out$clus <- apply(X = U, MARGIN = 1, FUN = which.max) # clusters
  return(out)
}

FKM.mahal.opt <- function(X, k, m = 2, stand = 0, conv = 1e-9, 
                          maxit = 1e6, prior, ntau = 1e3, 
                          eps = 1e-3, alpha){
  library(clue)
  # matrix of the true partition
  n <- nrow(X)
  q <- ncol(X)
  Ut <- matrix(0, n, k)
  for(i in 1:n){
    Ut[i, prior[i]] <- 1
  }
  # step 1
  d.star <- Inf
  ds.star <- rep(Inf, k)
  M <- list()
  for(j in 1:k){
    M[[j]] <- diag(1, q)
  }
  tau <- 1
  G.star <- list()
  G <- list()
  delta <- numeric(k)
  ds <- numeric(k)
  out <- list()
  # iterate until convergence or maximum iterations
  while(tau <= ntau & d.star >= eps){
    if(tau == 1){
      startU <- NULL
    }else{
      startU <- fkm$U
    }
    # step 2
    fkm <- FKM.mahal(X = X, k = k, m = m, stand = stand, startU = startU, 
                     conv = conv, maxit = maxit, clus_var = M)
    # reorder cluster indices
    if(tau == 1){
      count <- table(fkm$clus, prior)
      match <- solve_LSAP(count, maximum = TRUE)
      fkm$clus <- match[fkm$clus]
      fkm$U[,match] <- fkm$U
    }
    # step 3
    for(j in 1:k){
      ds[j] <- mean(abs(Ut[,j] - fkm$U[,j]))
    }
    d <- mean(ds)
    # step 4
    if(d < d.star){
      d.star <- d
      ds.star <- ds
      M.star <- M
    }
    if(d.star < eps || tau == ntau){
      out$d <- d.star
      out$ds <- ds.star
      out$M <- M.star
      out$FKM <- fkm
      break
    }
    # step 5
    tau <- tau + 1
    for(j in 1:k){
      G.star[[j]] <- chol(M.star[[j]]) # Cholesky decomposition
      G[[j]] <- G.star[[j]]
      delta[j] <- alpha * ds.star[j]
      tmp <- G.star[[j]] + matrix(rnorm(q*q, 0, delta[j]), q, q)
      G[[j]][upper.tri(G[[j]], diag = TRUE)] <- 
        tmp[upper.tri(tmp, diag = TRUE)]
      rm(tmp)
      diag(G[[j]])[diag(G[[j]]) < 10e-5] <- 10e-5
      M[[j]] <- t(G[[j]]) %*% G[[j]]
    }
  } 
  return(out)
}

dist_cluster <- function(prior, U){
  library(clue)
  k <- length(table(prior))
  n <- length(prior)
  Ut <- matrix(0, n, k)
  for(i in 1:n){
    Ut[i, prior[i]] <- 1
  }
  clust <- apply(X = U, MARGIN = 1, FUN = which.max)
  count <- table(clust, prior)
  match <- solve_LSAP(count, maximum = TRUE)
  clust <- match[clust]
  U[,match] <- U
  ds <- numeric(k)
  for(j in 1:k){
    ds[j] <- mean(abs(Ut[,j] - U[,j]))
  }
  d <- mean(ds)
  out <- list()
  out$d <- d
  out$ds <- ds
  return(out)  
}

