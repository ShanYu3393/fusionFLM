plot.beta <- function(fitted, file.name) {
  beta.path <- fitted$beta.path
  n.path <- nrow(beta.path)
  n.grid <- ncol(beta.path)
  subgroup <- rep(1:n.path, each = n.grid)
  y <- c(beta.path)
  x <- rep(1:n.grid, each = n.path)
  data.beta <- data.frame(x = x, y = y, subgroup = subgroup)
  p <- ggplot(data = data.beta, aes(x, y)) + 
    geom_line() + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  pdf(file = paste0("plots/", file.name, "beta_Lasso.pdf"), width = 3, height = 3, compress = TRUE)
  print(p)
  dev.off()
}
