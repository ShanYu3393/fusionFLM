admm.flm <- function(Z.list, Y, theta.init = NULL, eta.init = NULL, nu.init = NULL,
                 rho, lambda, nbasis = nbasis, edge.matrix, 
                 objective.path = FALSE, tau = 3.7, penalty = "Lasso") {
  gap <- c(1, rep(nbasis, length(Z.list) - 1))
  
  n.X <- length(Z.list) - 1
  n.group <- ncol(Z.list[[1]])
  if(is.matrix(edge.matrix)) {
    tmp <- list()
    for(iter in 1:(1+n.X)) tmp[[iter]] <- edge.matrix
    edge.matrix <- tmp   
  }
  
  n.group.diff <- sapply(edge.matrix, nrow)
  group <- c(1:n.group.diff[1], rep(1:sum(n.group.diff[-1]), each = nbasis) +
               n.group.diff[1]) 
  # notice that in this version, Z is a list.
  Z <- NULL
  nZ.all <- 0
  nA.all <- 0
  AA <- A <- list()
  
  iter.s <- 1
  iter.e <- length(gap)
  
  for (iter in iter.s:iter.e) {
    aa <- iter
    nZ <- ncol(Z.list[[aa]])
    nZ.all <- c(nZ.all, nZ)
    A[[iter]] <- kronecker(edge.matrix[[iter]], diag(gap[aa]))
    nA.all <- c(nA.all, dim(A[[iter]])[1])
    AA[[iter]] <- crossprod(A[[iter]])
    Z <- cbind(Z, Z.list[[aa]])
  }
  
  A <- bdiag(A)
  
  theta <- theta.init
  eta <- eta.init
  nu <- nu.init
  theta.norm <- NULL
  snorm <- 1
  resnorm <- 1
  if (is.null(theta.init)) theta <- rep(0, ncol(Z))
  if (is.null(eta.init)) eta <- rep(0, dim(A)[1])
  if (is.null(nu.init)) nu <- rep(0, dim(A)[1])
  
  res.norm <- 1
  if (length(lambda) == 1) {
    lambda <- rep(lambda, max(group))
  } else {
    if (length(lambda) == length(Z.list)) lambda <- rep(lambda, each = n.group.diff)
  }
  
  
  ZA.inver <- solve(crossprod(Z) + rho * bdiag(AA))
  YZ <- crossprod(Y, Z)
  Step <- 0
  # theta=1
  objective <- c(10, 100)
  theta.diff <- thetanorm <- 10
  while (Step <= 2000 & theta.diff > 0.01 * thetanorm) {
    Step <- Step + 1
    
    thetanorm <- sqrt(mean(theta^2))
    theta.prior <- theta
    
    # Step 2: update eta, Soft Thresholding
    W0 <- c()
    for (iter in 1:length(gap)) {
      a <- sum(nZ.all[1:iter]) + 1
      b <- sum(nZ.all[1:(iter + 1)])
      theta.matrix <- matrix(theta[a:b], nrow = gap[iter])
      W0.iter <- c(apply(edge.matrix[[iter]], 1, function(x) {
        theta.matrix[, which(x == 1)] - theta.matrix[, which(x == -1)]
      }))
      W0 <- c(W0, W0.iter)
    }
    W <- W0 + 1 / rho * nu
    
    
    W.norm <- tapply(W, group, function(x) sqrt(sum(x^2)))
    if (penalty == "Lasso") {
      a <- apply(
        cbind(W.norm, lambda), 1,
        function(x) drop(soft.thresh(x[1], x[2] / rho) / x[1])
      )
    }
    
    if (penalty == "SCAD") {
      a <- apply(
        cbind(W.norm, lambda), 1,
        function(x) drop(ST_SCAD(x[1], tau, rho, x[2]) / x[1])
      )
    }
    
    a[which(is.nan(a))] <- 0
    eta <- rep(a, times = table(group)) * W
    
    
    # Step 1: update theta, Quadratic Optimization
    theta <- ZA.inver %*% t(YZ - crossprod(nu - rho * eta, A))
    theta.norm <- c(theta.norm, as.matrix(sqrt(crossprod(theta))))
    res <- W0 - eta
    nu <- nu + rho * (res)
    theta.diff <- sum((theta.prior - theta)^2)
    
    # calculate objective value
    if (objective.path == TRUE) {
      obj.iter <- 0.5 * sum((Y - Z %*% theta)^2) + sum(lambda * tapply(A %*% theta, group, function(x) sqrt(crossprod(x))))
      objective <- c(objective, obj.iter)
    }
  }
  
  # loss
  Y.hat <- Z %*% theta
  loss <- sum((Y - Y.hat)^2)
  
  list(
    theta = theta, loss = loss, converge.step = Step, eta = eta, Yhat = Y.hat,
    objective = objective, theta.norm = theta.norm, resnorm = resnorm, snorm = snorm
  )
}
