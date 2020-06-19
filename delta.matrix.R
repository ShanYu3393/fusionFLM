delta.matrix <- function(n, gap = 1) {
  
  delta.matrix <- matrix(ncol = n, nrow = n * (n - 1) / 2)
  star <- 1
  end <- n - 1

  for (i in 1:(n - 1)) {
    if (i == 1) {
      delta.matrix[(star:end), ] <- cbind(1, -diag(n - i))
    } else {
      delta.matrix[(star:end), ] <- cbind(matrix(rep(0, (i - 1) * (n - i)),
                                                 nrow = n - i), rep(1, n - i), -diag(n - i))
    }

    star <- star + (n - i)
    end <- star + (n - i - 2)
  }

  kronecker(delta.matrix, diag(gap))
}
