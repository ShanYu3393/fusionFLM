kmeans.fusionFLM <- function (data, sp.basis) {
  
  nbasis <- sp.basis$nbasis
  n.X <- length(data$X.list)
  n <- length(data$Y)
  t.ober <- data$t.ober
  X.list <- data$X.list
  t.range <- sp.basis$rangeval
  Int <- data$subgroup.id
  Y <- data$Y
  n.group <- ncol(data$subgroup.id)
  
  # calculate W matrix ----------------------
  W <- matrix(ncol = (nbasis) * n.X, nrow = n)
  for (iter in 1:n.X) {
    W1 <- sapply(1:n, function(x) int.2(sp.basis, t.ober[[x]],
                                        X.list[[iter]][[x]], t.range))
    W[, ((iter - 1) * (nbasis) + 1):(iter * (nbasis))] <- t(W1)
  }
  
  # calculate Z matrix --------------------
  Z.list <- list()
  Z.list[[1]] <- Int
  n.group <- ncol(Int)
  for (iter in 1:n.X) {
    Z.iter <- matrix(0, ncol = (nbasis) * n.group, nrow = n)
    for (i in 1:n.group) {
      Z.iter[which(Int[, i] == 1), ((i - 1) * (nbasis) + 1):(i * (nbasis))] <-
        W[which(Int[, i] == 1), ((iter - 1) * (nbasis) + 1):(iter * (nbasis))]
    }
    Z.list[[(iter + 1)]] <- Z.iter
  }
  Z <- do.call('cbind', Z.list)
  
  # obtain individual estimators -----------
  theta.ind <- solve(crossprod(Z)) %*% t(Z) %*% Y
  beta.ind <- matrix(theta.ind[1:n.group], ncol = 1)
  gamma.ind <- matrix(theta.ind[-(1:n.group)], ncol = nbasis, byrow = TRUE)
  
  kmeans.beta.all <- clusGap(beta.ind, FUN = kmeans, K.max = 10)$Tab[, 3]
  kmeans.beta <- kmeans(beta.ind, which.max(kmeans.beta.all))
  group <- kmeans.beta$cluster
  theta.hat <- kmeans.beta$centers[kmeans.beta$cluster]
  for (iX in 1:n.X) {
    gamma.ind.iX <- gamma.ind[((iX - 1) * n.group + 1):(iX * n.group), ]
    kmeans.alpha.all <- clusGap(gamma.ind.iX, FUN = kmeans, K.max = 10)$Tab[, 3]
    kmeans.alpha <- kmeans(gamma.ind.iX, which.max(kmeans.alpha.all))
    group <- rbind(group, kmeans.alpha$cluster)
  }
  
  # create Q.list based on the estimated K-means group -----------------
  Q.list.Kmeans <- list()
  for(iter in 1:nrow(group)) {
    if(length(unique(group[iter, ])) == 1){
      Q.list.Kmeans[[iter]] <- matrix(group[iter, ], ncol = 1)
    } else {
      Q.list.Kmeans[[iter]] <- model.matrix(~as.factor(group[iter, ]) + 0)
    }
    
  }
  
  theta.hat <- oracle.FLM(data, Q.list = Q.list.Kmeans, sp.basis = sp.basis)$theta.hat
  
  list(group = group, theta.hat = theta.hat) 
}