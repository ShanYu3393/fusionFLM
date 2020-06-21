rm(list = ls())
# library
library(igraph)
library(fda)
library(fossil)
library(mclust)
library(cluster)
library(parallel)

# source files
source('alpha.func.5.R')
source('int.R')
source('int.2.R')
source('simu.data.generating.R')
source('delta.matrix.R')
source('edge.generating.R')
source('fit.fusionFLM.R')
source('plot.alpha.R')
source('plot.beta.R')
source('admm.flm.R')
source('group.identify.R')
source('fusionFLM.R')
source('soft.thresh.R')
source('jaccard.R')
source('eval.FLM.R')
source('kmeans.fusionFLM.R')
source('ind.FLM.R')
source('oracle.FLM.R')


# read in data
load('data/kinship.matrix.yield.rda')
load('data/true.group.5.rda')

# simulation setting
n.subgroup <- length(true.group.5)
n.ober <- 50
error.sd <- 0.5
number.knn <- 10
n.T <- 100
n.X <- 1
t.range <- c(0, 1)
Q.alpha1 <- model.matrix(~ as.factor(true.group.5) + 0)
Q.intercept <- cbind(Q.alpha1[, 1] + Q.alpha1[, 5], rowSums(Q.alpha1[, 2:4]))
Q.list <- list(Q.intercept, Q.alpha1)
beta <- c(5, 6)
alpha <- list(list(alpha1, alpha2, alpha3, alpha4, alpha5))
true.cluster = rbind(Q.intercept[, 1] + 2 * Q.intercept[, 2], true.group.5)
true.beta <- beta
true.alpha <- alpha

# true coefficient function
col <- categorical_pal(5)
time.point <- (0:100) * 0.01
plot(time.point, alpha1(time.point), type = "l", col = col[1], lwd = 2,
     ylab = expression(alpha), xlab = 't')
points(time.point, alpha2(time.point), type = "l", col = col[2], lwd = 2)
points(time.point, alpha3(time.point), type = "l", col = col[3], lwd = 2)
points(time.point, alpha4(time.point), type = "l", col = col[4], lwd = 2)
points(time.point, alpha5(time.point), type = "l", col = col[5], lwd = 2)

# spline basis
N <- 4
rho <- 2
nbasis <- N + rho + 1
sp.basis <- create.bspline.basis(t.range, nbasis = nbasis, norder = rho + 1)

# construct graph for the fused lasso
edge.matrix.prior <- edge.generating(subgroup.name = 1:n.subgroup,
                                     type = 'knn', dist = kinship.matrix.yield,
                                     knn = number.knn, vertex.size = 6,
                                     Group.est = true.group.5)

# fit the model
lambda1 <- 10^seq(-2, 2, by = 1)
lambda2 <- 10^seq(-2, 2, by = 1)
Lambda1 <- expand.grid(lambda1, lambda2)

lambda1 <- 10^seq(-2, 2, by = 1)
lambda2 <- 10^seq(-2, 2, by = 1)
Lambda2 <- expand.grid(lambda1, lambda2)

Lambda.list = list(Lambda1, Lambda2)

# generate data 
simu.main <- function (iter) {
  
  set.seed(iter)
  print(iter)
  data.simu <- simu.data.generating(n.subgroup, n.ober, n.T, n.X,
                                    Q.list, beta, alpha, error.sd, t.range = t.range)
  
  t0 <- proc.time()
  fitted.oracle <- oracle.FLM(data = data.simu, Q.list = Q.list,
                              sp.basis = sp.basis)
  t.oracle <- proc.time() - t0 
  
  # individual group estimator
  t0 <- proc.time()
  fitted.ind <- ind.FLM(data = data.simu, sp.basis)
  t.ind <- proc.time() - t0  
  
  # functional linear regression, k-means
  t0 <- proc.time()
  fitted.kmeans <- kmeans.fusionFLM(data = data.simu, sp.basis)
  t.kmeans <- proc.time() - t0
  
  # functional linear regression, sparse graph
  t0 <- proc.time()
  fitted.sparse.g1 <- fusionFLM(data = data.simu, sp.basis = sp.basis, 
                                edge.matrix = edge.matrix.prior,
                                Lambda.list = Lambda.list, initial.type = 'kmeans', 
                                objective.path = FALSE, save.plot = FALSE)
  t.sparse.g1 <- proc.time() - t0
  
  t0 <- proc.time()
  fitted.sparse.g2 <- fusionFLM(data = data.simu, sp.basis = sp.basis, 
                                edge.matrix = edge.matrix.prior,
                                Lambda.list = Lambda.list, initial.type = 'individual', 
                                objective.path = FALSE, save.plot = FALSE)
  t.sparse.g2 <- proc.time() - t0
  
  t0 <- proc.time()
  fitted.sparse.g3 <- fusionFLM(data = data.simu,  sp.basis = sp.basis, 
                                edge.matrix = edge.matrix.prior,
                                Lambda.list = Lambda.list, initial.type = 'lasso', 
                                objective.path = FALSE, save.plot = TRUE)
  t.sparse.g3 <- proc.time() - t0
  
  result.oracle <- eval.FLM(fitted.oracle, true.cluster, true.beta, true.alpha, sp.basis)
  result.ind <- eval.FLM(fitted.ind, true.cluster, true.beta, true.alpha, sp.basis)
  result.kmeans <- eval.FLM(fitted.kmeans, true.cluster, true.beta, true.alpha, sp.basis)
  result.sparse.g1 <- eval.FLM(fitted.sparse.g1, true.cluster, true.beta, true.alpha, sp.basis)
  result.sparse.g2 <- eval.FLM(fitted.sparse.g2, true.cluster, true.beta, true.alpha, sp.basis)
  result.sparse.g3 <- eval.FLM(fitted.sparse.g3, true.cluster, true.beta, true.alpha, sp.basis)
  
  time.all <- c(t.oracle[3], t.ind[3], t.kmeans[3], t.sparse.g1[3], t.sparse.g2[3],
                t.sparse.g3[3])
  cluster.all <- rbind(result.oracle$cluster.results,
                       result.ind$cluster.results,
                       result.kmeans$cluster.results,
                       result.sparse.g1$cluster.results,
                       result.sparse.g2$cluster.results,
                       result.sparse.g3$cluster.results)
  rmse.beta.all <- rbind(result.oracle$rmse.beta,
                         result.ind$rmse.beta,
                         result.kmeans$rmse.beta,
                         result.sparse.g1$rmse.beta,
                         result.sparse.g2$rmse.beta,
                         result.sparse.g3$rmse.beta)
  rmise.alpha.all <- rbind(result.oracle$rmise.alpha,
                           result.ind$rmise.alpha,
                           result.kmeans$rmise.alpha,
                           result.sparse.g1$rmise.alpha,
                           result.sparse.g2$rmise.alpha,
                           result.sparse.g3$rmise.alpha)
  list(time.all, cluster.all, rmse.beta.all, rmise.alpha.all)
}

# result <- simu.main(1)
# print(result)
result <- mclapply(1:100, simu.main, mc.cores = 16)
save(file = paste0('result/cluster5', 'n', n.ober, 'sigma', error.sd, '.rda'), result)