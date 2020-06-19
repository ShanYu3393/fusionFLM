# given group structure and obtain estimators
oracle.FLM <- function (data, Q.list, sp.basis) {
  
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
  Z.list.oracle <- list()
  Z.list.oracle[[1]] <- Int %*% Q.list[[1]]
  n.group <- ncol(Int)
  for (iter in 1:n.X) {
    Z.iter <- matrix(0, ncol = (nbasis) * n.group, nrow = n)
    for (i in 1:n.group) {
      Z.iter[which(Int[, i] == 1), ((i - 1) * (nbasis) + 1):(i * (nbasis))] <-
        W[which(Int[, i] == 1), ((iter - 1) * (nbasis) + 1):(iter * (nbasis))]
    }
    Z.list.oracle[[(iter + 1)]] <- Z.iter %*% kronecker(Q.list[[iter+1]], diag(nbasis))
  }
  Z.oracle <- do.call('cbind', Z.list.oracle)
  
  # coefficient estimation --------------------------------------
  theta.oracle <- solve(crossprod(Z.oracle)) %*% t(Z.oracle) %*% Y
  
  # expand to full coefficients
  theta.hat <- Q.list[[1]] %*% matrix(theta.oracle[1:ncol(Q.list[[1]])], ncol = 1)
  gamma.oracle <- matrix(theta.oracle[-(1:ncol(Q.list[[1]]))], ncol = nbasis, byrow = TRUE)
  
  start <- 1
  end <- ncol(Q.list[[2]])
  for (iX in 1:n.X) {
    gamma.iX <- t(Q.list[[iX+1]] %*% gamma.oracle[start:end, ])
    theta.hat <- c(theta.hat, c(gamma.iX))
    if(n.X > 1) {
      start <- ncol(Q.list[[iX+1]]) + start
      end <- end + ncol(Q.list[[iX+2]])
    }
  }
 
  group.list <- lapply(Q.list, function (x) 
    apply(x, 1, function (tmp) which(tmp == 1)))
  group <- do.call('rbind', group.list)
  
  list(group = group, theta.hat = theta.hat) 
}
