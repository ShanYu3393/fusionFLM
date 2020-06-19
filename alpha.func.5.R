# true coefficient functions
alpha1 <- function(x) cos(2 * pi * x)
alpha1 <- Vectorize(alpha1)

alpha2 <- function(x) sin(2 * pi * x)
alpha2 <- Vectorize(alpha2)

alpha3 <- function(x) 1.5 * (x - 0.5)
alpha3 <- Vectorize(alpha3)

alpha4 <- function(x) 1 - 2 * exp(-6 * x)
alpha4 <- Vectorize(alpha4)

alpha5 <- function(x) 2 * exp(-6 * x) - 1
alpha5 <- Vectorize(alpha5)


