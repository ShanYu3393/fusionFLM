fusionFLM <- function (data, sp.basis = NULL,
                       sp.basis1 = NULL, sp.basis2 = NULL, 
                       edge.matrix, Lambda.list,
                       initial.type = 'lasso', 
                       objective.path = FALSE, save.plot) {

  if (!is.null(sp.basis) & is.null(sp.basis1) & is.null(sp.basis2)) {
    sp.basis1 <- sp.basis
    sp.basis2 <- sp.basis
  } else if (is.null(sp.basis) & !is.null(sp.basis1) & !is.null(sp.basis2)) {
    sp.basis1 <- sp.basis1
    sp.basis2 <- sp.basis1
  } else {
    stop('no valid spline basis functions.')
  }
  
  # calculate some basic parameters ---------------
  n.X <- length(data$X.list)
  n <- length(data$Y)
  t.ober <- data$t.ober
  X.list <- data$X.list
  t.range <- sp.basis1$rangeval
  Int <- data$subgroup.id
  Y <- data$Y
  n.group <- ncol(data$subgroup.id)
  
  if(is.matrix(edge.matrix)) {
    tmp <- list()
    for(iter in 1:(1+n.X)) tmp[[iter]] <- edge.matrix
    edge.matrix <- tmp   
  }
  n.group.diff <- sapply(edge.matrix, nrow)
  
  
  
  # initial estimators -------
  nbasis1 <- sp.basis1$nbasis
  
  # calculate W matrix 
  W.init <- matrix(ncol = (nbasis1) * n.X, nrow = n)
  for (iter in 1:n.X) {
    W1 <- sapply(1:n, function(x) int.2(sp.basis1, t.ober[[x]],
                                        X.list[[iter]][[x]], t.range))
    W.init[, ((iter - 1) * (nbasis1) + 1):(iter * (nbasis1))] <- t(W1)
  }
  
  # calculate Z matrix 
  Z.list.init <- list()
  Z.list.init[[1]] <- Int
  n.group <- ncol(Int)
  for (iter in 1:n.X) {
    Z.iter <- matrix(0, ncol = (nbasis1) * n.group, nrow = n)
    for (i in 1:n.group) {
      Z.iter[which(Int[, i] == 1), ((i - 1) * (nbasis1) + 1):(i * (nbasis1))] <-
        W.init[which(Int[, i] == 1), ((iter - 1) * (nbasis1) + 1):(iter * (nbasis1))]
    }
    Z.list.init[[(iter + 1)]] <- Z.iter
  }
  Z.init <- do.call('cbind', Z.list.init)
  
  # obtain individual estimators 
  theta.ind.init <- solve(crossprod(Z.init)) %*% t(Z.init) %*% Y
  gamma.ind.init <- matrix(theta.ind.init[-(1:n.group)], ncol = nbasis1, byrow = TRUE)
  
  # group structure based on the initial stage.
  group.init <- c(1:n.group.diff[1], rep(1:sum(n.group.diff[-1]), each = nbasis1) +
               n.group.diff[1])
  
  # adaptive lasso estimators -------
  nbasis2 <- sp.basis2$nbasis
  
  # calculate W matrix 
  W.adp <- matrix(ncol = (nbasis2) * n.X, nrow = n)
  for (iter in 1:n.X) {
    W1 <- sapply(1:n, function(x) int.2(sp.basis2, t.ober[[x]],
                                        X.list[[iter]][[x]], t.range))
    W.adp[, ((iter - 1) * (nbasis2) + 1):(iter * (nbasis2))] <- t(W1)
  }
  
  # calculate Z matrix 
  Z.list.adp <- list()
  Z.list.adp[[1]] <- Int
  n.group <- ncol(Int)
  for (iter in 1:n.X) {
    Z.iter <- matrix(0, ncol = (nbasis2) * n.group, nrow = n)
    for (i in 1:n.group) {
      Z.iter[which(Int[, i] == 1), ((i - 1) * (nbasis2) + 1):(i * (nbasis2))] <-
        W.adp[which(Int[, i] == 1), ((iter - 1) * (nbasis2) + 1):(iter * (nbasis2))]
    }
    Z.list.adp[[(iter + 1)]] <- Z.iter
  }
  Z.adp <- do.call('cbind', Z.list.adp)
  
  # obtain individual estimators 
  theta.ind.adp <- solve(crossprod(Z.adp)) %*% t(Z.adp) %*% Y
  gamma.ind.adp <- matrix(theta.ind.adp[-(1:n.group)], ncol = nbasis2, byrow = TRUE)
  
  # group structure based on the initial stage.
  group.adp <- c(1:n.group.diff[1], rep(1:sum(n.group.diff[-1]), each = nbasis2) +
                    n.group.diff[1])
  
  # calculate estimators -------------------------------------
  if(initial.type == 'lasso'){
    fitted <- fit.fusionFLM(Z.list.init, Y, Lambda = Lambda.list[[1]], 
                            edge.matrix, theta.init = theta.ind.init, weight = 1,
                            rho = 1, nbasis = nbasis1,
                            objective.path = FALSE, criterion = 'BIC',
                            save.plot = save.plot,
                            file.name = 'alpha.lasso')
    
    if(save.plot == TRUE) {
      plot.beta(fitted, file.name = 'beta.lasso')
    }
    
    weight <- 1 / tapply(fitted$eta.hat, group.init, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list.adp, Y, Lambda = Lambda.list[[2]],
                            edge.matrix, theta.init = theta.ind.adp,
                            weight = weight,
                            rho = 1, nbasis = nbasis2,
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
      edge.matrix[[1]] %*% theta.ind.init[1:n.group],
      as.vector(bdiag(lapply(edge.matrix, 
                   function(x) kronecker(x, diag(nbasis1)))[[-1]]) %*%
        theta.ind.init[-(1:n.group)])
    )
    
    weight <- 1 / tapply(eta.init, group.init, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list.adp, Y, Lambda = Lambda.list[[2]], 
                            edge.matrix, theta.init = theta.ind.adp, weight = weight,
                            rho = 1, nbasis = nbasis2,
                            objective.path = objective.path,
                            criterion = 'EBIC',
                            save.plot = save.plot)
  }
  
  if (initial.type == 'kmeans'){
    fitted.kmeans <- kmeans.fusionFLM(data = data, sp.basis1)
    theta.init <- fitted.kmeans$theta.hat
    eta.init <- c(
      edge.matrix[[1]] %*% theta.init[1:n.group],
      as.vector(bdiag(lapply(edge.matrix, 
                   function(x) kronecker(x, diag(nbasis1)))[[-1]]) %*%
        theta.init[-(1:n.group)])
    )
    
    weight <- 1 / tapply(eta.init, group.init, function(x) sqrt(crossprod(x)))
    weight <- weight^2
    weight[which(weight > 1e+06)] <- 1e+06
    
    fitted <- fit.fusionFLM(Z.list.adp, Y, Lambda = Lambda.list[[2]], 
                            edge.matrix, theta.init = theta.ind.adp, weight = weight,
                            rho = 1, nbasis = nbasis2,
                            objective.path = objective.path,
                            criterion = 'EBIC',
                            save.plot = save.plot)
  }
  fitted
}