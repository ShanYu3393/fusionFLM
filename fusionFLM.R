fusionFLM <- function (data, sp.basis, edge.matrix, Lambda.list,
                       initial.type = 'lasso', 
                       objective.path = FALSE, save.plot) {
  
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
  gamma.ind <- matrix(theta.ind[-(1:n.group)], ncol = nbasis, byrow = TRUE)
  
  if(is.matrix(edge.matrix)) {
    tmp <- list()
    for(iter in 1:(1+n.X)) tmp[[iter]] <- edge.matrix
    edge.matrix <- tmp   
  }
  
  n.group.diff <- sapply(edge.matrix, nrow)
  group <- c(1:n.group.diff[1], rep(1:sum(n.group.diff[-1]), each = nbasis) +
               n.group.diff[1])
  
  if(initial.type == 'lasso'){
    fitted <- fit.fusionFLM(Z.list, Y, Lambda = Lambda.list[[1]], 
                            edge.matrix, theta.init = theta.ind, weight = 1,
                            rho = 1, nbasis = nbasis,
                            objective.path = FALSE, criterion = 'BIC',
                            save.plot = save.plot,
                            file.name = 'alpha.lasso')
    
    if(save.plot == TRUE) {
      plot.beta(fitted, file.name = 'beta.lasso')
    }
    
    weight <- 1 / tapply(fitted$eta.hat, group, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list, Y, Lambda = Lambda.list[[2]],
                            edge.matrix, theta.init = theta.ind,
                            weight = weight,
                            rho = 1, nbasis = nbasis,
                            objective.path = objective.path,
                            criterion = 'EBIC',
                            save.plot = save.plot,
                            file.name = 'alpha.alasso')
    
    if(save.plot == TRUE) {
      plot.beta(fitted, file.name = 'beta.alasso')
    }
  }
  
  if(initial.type == 'individual'){
    
    eta.init <- c(
      edge.matrix[[1]] %*% theta.ind[1:n.group],
      as.vector(bdiag(lapply(edge.matrix, 
                   function(x) kronecker(x, diag(nbasis)))[[-1]]) %*%
        theta.ind[-(1:n.group)])
    )
    
    weight <- 1 / tapply(eta.init, group, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list, Y, Lambda = Lambda.list[[2]], 
                            edge.matrix, theta.init = theta.ind, weight = weight,
                            rho = 1, nbasis = nbasis,
                            objective.path = objective.path,
                            criterion = 'EBIC',
                            save.plot = save.plot)
  }
  
  if (initial.type == 'kmeans'){
    beta.ind <- matrix(theta.ind[1:n.group], ncol = 1)
    gamma.ind <- matrix(theta.ind[-(1:n.group)], ncol = nbasis, byrow = TRUE)
    
    kmeans.beta.all <- clusGap(beta.ind, FUN = kmeans, K.max = 10)$Tab[, 3]
    kmeans.beta <- kmeans(beta.ind, which.max(kmeans.beta.all))
    theta.init <- kmeans.beta$centers[kmeans.beta$cluster]
    for (iX in 1:n.X) {
      gamma.ind.iX <- gamma.ind[((iX - 1) * n.group + 1):(iX * n.group), ]
      kmeans.alpha.all <- clusGap(gamma.ind.iX, FUN = kmeans, K.max = 10)$Tab[, 3]
      kmeans.alpha <- kmeans(gamma.ind.iX, which.max(kmeans.alpha.all))
      theta.init <- c(theta.init, c(t(kmeans.alpha$centers[kmeans.alpha$cluster, ])))
    }
    eta.init <- c(
      edge.matrix[[1]] %*% theta.init[1:n.group],
      as.vector(bdiag(lapply(edge.matrix, 
                   function(x) kronecker(x, diag(nbasis)))[[-1]]) %*%
        theta.init[-(1:n.group)])
    )
    
    weight <- 1 / tapply(eta.init, group, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list, Y, Lambda = Lambda.list[[2]], 
                            edge.matrix, theta.init = theta.ind, weight = weight,
                            rho = 1, nbasis = nbasis,
                            objective.path = objective.path,
                            criterion = 'EBIC',
                            save.plot = save.plot)
  }
  fitted
}