rm(list = ls())
# library
library(igraph)
library(fda)
library(fossil)
library(mclust)
library(cluster)
library(parallel)

# source files
source("alpha.func.5.R")
source("int.R")
source("int.2.R")
source("simu.data.generating.R")
source("delta.matrix.R")
source("edge.generating.R")
source("fit.fusionFLM.R")
source("plot.alpha.R")
source("admm.flm.R")
source("group.identify.R")
source("fusionFLM.R")
source("soft.thresh.R")
source("jaccard.R")
source("eval.FLM.R")
source("kmeans.fusionFLM.R")
source("ind.FLM.R")


# read in data
load("data/kinship.sub.yield.rda")
load("data/yield.sub.rda")
load("data/covariate.sub.yield.rda")


# spline basis
N <- 2
rho <- 2
nbasis <- N + rho + 1
t.range <- c(0, 2)
sp.basis <- create.bspline.basis(t.range, nbasis = nbasis, norder = rho + 1)

n.subgroup <- ncol(kinship.sub)

# construct graph for the fused lasso
edge.matrix.prior <- edge.generating(
  subgroup.name = 1:n.subgroup,
  type = "mst", dist = kinship.sub
)

# edge.matrix.prior <- edge.generating(
#   subgroup.name = 1:n.subgroup,
#   type = "knn", dist = kinship.sub,
#   knn = 8, vertex.size = 6,
#   Group.est = 1
# )
# 
# edge.matrix.prior <- edge.generating(
#   subgroup.name = 1:n.subgroup,
#   type = "threshhold", dist = kinship.sub, cutoff.value = 0.2
# )

# fit the model
# fused lasso estimation
lambda1 <- 10^seq(0, 3, by = 1)
lambda2 <- 10^seq(0, 2, by = 0.5)
Lambda1 <- expand.grid(lambda1, lambda2)
# Lambda1=c(100,100)
lambda1 <- 10^seq(-2, 2, by = 1)
lambda2 <- 10^seq(-2, 0, by = 1)
Lambda2 <- expand.grid(lambda1, lambda2)

Lambda.list <- list(Lambda1, Lambda2)

# generate data
dat.yield <- list()
dat.yield$Y <- yield.sub$Yield
dat.yield$X.list <- list(covariate_ts$SR)
dat.yield$t.ober <- covariate_ts$Time
dat.yield$subgroup.id <- as.matrix(model.matrix(~ yield.sub$Pedigree + 0))

# # individual group estimator
t0 <- proc.time()
fitted.ind <- ind.FLM(data = dat.yield, sp.basis)
t.ind <- proc.time() - t0
# 
# # functional linear regression, k-means
# t0 <- proc.time()
# fitted.kmeans <- kmeans.fusionFLM(data = dat.yield, sp.basis)
# t.kmeans <- proc.time() - t0

# functional linear regression, sparse graph
t0 <- proc.time()
fitted.sparse <- fusionFLM(
  data = dat.yield, sp.basis, edge.matrix.prior,
  Lambda.list, initial.type = "lasso",
  objective.path = FALSE, save.plot = FALSE
)
t.sparse <- proc.time() - t0
fitted.sparse <- fitted.sparse.g1
save(file = 'result/fitted.application.7.RData', fitted.sparse)

