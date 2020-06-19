int <- function(func, time, X, t.range) {
  sum(X * func(time) * (c(time[-1], t.range[2]) - time))
}
