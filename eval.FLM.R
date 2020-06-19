eval.FLM <- function(fitted, true.cluster, true.beta, true.alpha,
                           sp.basis) {
  
  # clustering results
  cluster.results <- c()
  for (iter in 1:nrow(true.cluster)) {
    RI <- rand.index(true.cluster[iter, ], fitted$group[iter, ])
    aRI <- adjustedRandIndex(fitted$group[iter, ], true.cluster[iter, ])
    Jcard <- jaccard(true.cluster[iter, ], fitted$group[iter, ])
    group.size <- max(fitted$group[iter, ])
    cluster.results <- c(cluster.results, c(RI, aRI, Jcard, group.size))
  }
    
    
  # RMSE of beta
  n.group <- ncol(true.cluster)
  if(length(true.beta) == 1) {
    beta.all = true.beta
  } else {
    beta.all <- model.matrix(~ as.factor(true.cluster[1, ]) + 0) %*% 
      true.beta
  }
  
  mse.beta.all <- (fitted$theta.hat[(1:n.group)] - beta.all)^2
  rmse.beta <- sqrt(tapply(mse.beta.all, true.cluster[1, ], mean))
  
  # RMISE of alpha
  t.range <- sp.basis$rangeval
  time <- seq(t.range[1], t.range[2], by = 0.01)
  gamma.all.matrix <- matrix(fitted$theta.hat[-(1:n.group)], ncol = n.group)
  rmise.alpha <- c()
  for(iter.X in 1:(nrow(true.cluster) - 1)){
    gamma.all.matrix.iX <- 
      gamma.all.matrix[, ((iter.X - 1) * n.group + 1):(iter.X * n.group)]
    mise.alpha.iX <- c()
    for (iter in 1:n.group) {
      index.alpha <- true.cluster[iter.X + 1,iter]
      alpha.time.true <- true.alpha[[iter.X]][[index.alpha]](time)
      alpha.time.est <- eval.basis(time, basisobj = sp.basis) %*% 
        gamma.all.matrix.iX[, iter]
      mise.alpha.iX <- c(mise.alpha.iX, 
                          sum(0.01 * (alpha.time.true - alpha.time.est)^2))
    }
    rmise.alpha <- c(rmise.alpha,
                  sqrt(tapply(mise.alpha.iX, true.cluster[iter.X + 1, ], mean)))
  }
  
  list(cluster.results = cluster.results, rmse.beta = rmse.beta,
       rmise.alpha = rmise.alpha)
}
