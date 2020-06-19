plot.alpha <- function(Coeff, basis, 
                       group = rep(1, ncol(Coeff)), lwd = 0.5) {
  # decide ylim
  t.range <- sp.basis$rangeval
  time <- seq(t.range[1], t.range[2], by = 0.01)
  value.alpha <- c(eval.basis(time, basisobj = sp.basis) %*% Coeff)
  ylim = c(min(value.alpha), max(value.alpha))
  
  # plot alpha functions
  col <- categorical_pal(max(group))
  plot(fd(coef = Coeff[, 1], basisobj = sp.basis),
       ylim = ylim, col = col[group[1]], lwd = lwd)
  for (iter in 2:ncol(Coeff)) {
    plot(fd(coef = Coeff[, iter], basisobj = sp.basis),
      add = TRUE, col = col[group[iter]],
      lwd = lwd
    )
  }
}
