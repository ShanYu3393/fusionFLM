soft.thresh <- function(x, tau) {
  if (tau <= 0) print("no valid tau")

  if (x > tau) {
    y <- x - tau
  } else if (x < -tau) {
    y <- x + tau
  } else {
    y <- 0
  }

  y
}
