#' Generating simulation examples
#' @para n.subgroup Number of subgroups K.
#' @para n.ober Number of observations within each subgroup.
#' @para n.T Number of observed time points.
#' @para t.range The endpoints of the time domain. Default is \code{c(0,1)}.
#' @para n.X Number of functional covariates.
#' @para Q.list A list contains matrices Q, which indicate the cluster
#' structures of subgroups.
#' @para beta A vector, which gives the true intercept values in each cluster.
#' @para alpha A list of functions, which gives the true coefficient functions
#' in each cluster.
#' @para error.sd Standard deviation of the noise.
#' @return A list of data information.
#' \item{Y}{response variable}
#' \item{X.list}{list of functional covariates.}
#' \item{t.ober}{observed time points for each subject.}
#' \item{subgroup.id}{a matrix indiates the subgroup of eahc subject.}
################################################################################
simu.data.generating <- function(n.subgroup, n.ober, n.T, n.X, Q.list, beta,
                                 alpha, error.sd, t.range = c(0, 1)) {
  
  subgroup.id = kronecker(diag(n.subgroup), matrix(rep(1, n.ober), ncol = 1))
  
  # generate random time points
  t.ober <- lapply(1:(n.subgroup * n.ober),
                   FUN = function(x) sort(runif(n.T)))

  # generate functional covariates X(t)
  X.list <- list()
  # Define the basis function of X(t)
  nbasis <- 9
  fbasis <- create.fourier.basis(t.range, nbasis)
  BX <- lapply(t.ober, FUN = eval.basis, basisobj = fbasis)
  
  for (iter in 1:n.X) {
    coeff <- rnorm(n.subgroup * n.ober * nbasis)
    coeff[coeff > 3] <- 3
    coeff[coeff < -3] <- -3
    coeff.matrix <- matrix(coeff, ncol = nbasis)
    X.list[[iter]] <- lapply(1:(n.subgroup * n.ober),
                             function(x) BX[[x]] %*% coeff.matrix[x, ])
  }

  # generate response variable within each cluster
  alpha.all <- list()
  a.all <- c()
  for (iter in 1:n.X) {
    # index of coefficient functions to use for each observation.
    index = apply(subgroup.id %*% Q.list[[iter+1]], 1,
                  function(x) which(x == 1))
    # integration of X(t)alpha(t) 
    a <- sapply(1:(n.subgroup * n.ober), function(x) 
      int(alpha[[iter]][[index[x]]], t.ober[[x]], X.list[[iter]][[x]], t.range))
    a.all <- cbind(a.all, a)
  }
  
  intercept <- subgroup.id %*% Q.list[[1]] %*% matrix(beta, ncol = 1)
  error <- rnorm(n.subgroup * n.ober, sd = error.sd)
  Y <- rowSums(a.all) + error + intercept - mean(error)

  list(Y = Y, X.list = X.list, t.ober = t.ober, subgroup.id = subgroup.id)
}
