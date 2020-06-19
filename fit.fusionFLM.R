fit.fusionFLM <- function(Z.list, Y, Lambda, edge.matrix,
                          theta.init, weight, rho = 1,
                          nbasis = nbasis, objective.path = FALSE,
                          save.plot = FALSE, criterion = 'EBIC',
                          file.name = NULL) {
  
  # some basis parameters
  n.group <- ncol(Z.list[[1]])
  n <- length(Y)
  n.X <- length(Z.list) - 1
  
  if(is.matrix(edge.matrix)) {
    tmp <- list()
    for(iter in 1:(1+n.X)) tmp[[iter]] <- edge.matrix
    edge.matrix <- tmp   
  }
  
  n.group.diff <- sapply(edge.matrix, nrow)
  group <- c(1:n.group.diff[1], rep(1:sum(n.group.diff[-1]), each = nbasis) +
               n.group.diff[1])
  Z.all <- do.call('cbind', Z.list)
  
  # initial value
  eta.init <- c(
    edge.matrix[[1]] %*% theta.init[1:n.group],
    as.vector(bdiag(lapply(edge.matrix, 
                 function(x) kronecker(x, diag(nbasis)))[[-1]]) %*%
      theta.init[-(1:n.group)])
  )
  
  theta.hat <- theta.init
  eta.hat <- eta.init
  nlambda <- nrow(Lambda)
  BIC <- rep(0, nlambda)
  THETA <- NULL
  beta.path <- c()
  CLUSTER <- list()
  
  for (iter in 1:(1 + n.X)) {
    CLUSTER[[iter]] <- matrix(0, nrow = nlambda, ncol = n.group)
  } 
  
  for (iter in 1:nlambda) {
    lambda <- c(rep(Lambda[iter, 1], n.group.diff[1]),
                rep(Lambda[iter, 2], sum(n.group.diff[-1]))) * weight
    
    cat(iter, "\n")
    result <- admm.flm(Z.list, Y,
                    theta.init = theta.hat, eta.init = eta.hat, nu.init = NULL,
                    rho, lambda, nbasis = nbasis, edge.matrix, 
                    objective.path = objective.path
    )
    coeff.hat <- list()
    coeff.hat[[1]] <- result$theta[1:n.group]
    for (jj in 2:(1 + n.X)) {
      coeff.hat[[jj]] <- result$theta[((n.group + 1) + (jj - 2) * nbasis * n.group):(n.group + (jj - 1) * nbasis * n.group)]
    }
    coeff.adj <- group.identify(result$eta, coeff.hat,
                             c(1, rep(nbasis, n.X)), edge.matrix)
    loss.new <- sum((Y - Z.all %*% coeff.adj$coeff.hat)^2)
    K <- sum(c(1, rep(nbasis, n.X)) * coeff.adj$cluster.num)
    eta.hat <- c(
      edge.matrix[[1]] %*% coeff.adj$coeff.hat[1:n.group],
      as.vector(bdiag(lapply(edge.matrix, 
                   function(x) kronecker(x, diag(nbasis)))[[-1]]) %*%
        coeff.adj$coeff.hat[-c(1:n.group)])
    )
    
    # calculate BIC
    if (criterion == 'BIC') {
      BIC[iter] <- log(loss.new / n) + log(n) * (K) / n
    } 
    if (criterion == 'EBIC') {
      BIC[iter] <- log(loss.new / n) + log(n) * (K) / n + log(length(result$theta)) * (K) / n
    }
    
    THETA <- cbind(THETA, coeff.adj$coeff.hat)
    for (ii in 1:(1 + n.X)) {
      CLUSTER[[ii]][iter, ] <- coeff.adj$cluster.grp[ii, ]
    }
    beta.path <- cbind(beta.path, coeff.adj$coeff.hat[(1:n.group)])
    
    if (save.plot == TRUE) {
      gamma.hat.matrix <- matrix(coeff.adj$coeff.hat[-(1:n.group)], nrow = nbasis)
      pdf(file = paste0("plots/", file.name, iter, ".pdf"), width = 3, height = 3, compress = TRUE)
      print(plot.alpha(gamma.hat.matrix, sp.basis))
      dev.off()
    }
  }
  plot(BIC)
  
  select.index <- which.min(BIC)
  theta.hat <- THETA[, select.index]
  selected.lambda <- Lambda[select.index, ]
  eta.hat <- c(
    edge.matrix[[1]] %*% as.vector(theta.hat[1:n.group]),
    as.vector(bdiag(lapply(edge.matrix, 
                 function(x) kronecker(x, diag(nbasis)))[[-1]]) %*%
      theta.hat[-c(1:n.group)])
  )
  
  group <- c()
  for (ii in 1:(1 + n.X)) {
    group <- rbind(group, CLUSTER[[ii]][select.index, ])
  }
  
  list(theta.hat = theta.hat, eta.hat = eta.hat, best.lam = selected.lambda, 
       group = group, beta.path = beta.path)
}
