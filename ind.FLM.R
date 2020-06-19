ind.FLM <- function (data, sp.basis) {
  
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
  theta.hat <- solve(crossprod(Z)) %*% t(Z) %*% Y
  group <- matrix(rep(1:n.group, each = 1 + n.X), nrow = n.X + 1)
  
  list(group = group, theta.hat = theta.hat) 
}