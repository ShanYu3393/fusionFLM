int.2 <- function(sp.basis, time, X, t.range) {
  time <- as.vector(time)
  X <- as.matrix(X)
  func.value <- t(eval.basis(time, basisobj = sp.basis))
  t(X * (c(time[-1], t.range[2]) - time)) %*% t(func.value)
}
