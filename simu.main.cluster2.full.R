rm(list = ls())

library(igraph)
library(fda)
library(fossil)
library(mclust)
library(cluster)
library(parallel)

# simulation generating
source('alpha.func.5.R')
source('simu.data.generating.R')
# basis functions in fusionFLM
source('int.R')
source('int.2.R')
source('delta.matrix.R')
source('edge.generating.R')
source('admm.flm.R')
source('group.identify.R')
source('soft.thresh.R')
# evaluation
source('jaccard.R')
source('eval.FLM.R')
source('plot.alpha.R')
# model fitting
source('fusionFLM.R')
source('fit.fusionFLM.R')
source('kmeans.fusionFLM.R')
source('ind.FLM.R')
source('oracle.FLM.R')

n.subgroup <- 40
n.ober <- 50
n.T <- 100
n.X <- 1
t.range <- c(0, 1)
Q.list <- list(matrix(rep(1, n.subgroup), ncol = 1),
               kronecker(diag(2), matrix(rep(1, n.subgroup/2), ncol = 1)))
beta <- 5
alpha <- list(list(alpha1, alpha2))
error.sd <- 0.5

# true coefficient function
col <- categorical_pal(5)
time.point <- (0:100) * 0.01
plot(time.point, alpha1(time.point), type = "l", col = col[1], lwd = 2,
     ylab = expression(alpha), xlab = 't')
points(time.point, alpha2(time.point), type = "l", col = col[2], lwd = 2)

# spline basis
N <- 4
rho <- 2
nbasis <- N + rho + 1
sp.basis <- create.bspline.basis(t.range, nbasis = nbasis, norder = rho + 1)

# fit the model
lambda1 <- 10^seq(-2, 2, by = 1)
lambda2 <- 10^seq(-2, 2, by = 1)
Lambda1 <- expand.grid(lambda1, lambda2)

lambda1 <- 10^seq(-2, 2, by = 1)
lambda2 <- 10^seq(-2, 2, by = 1)
Lambda2 <- expand.grid(lambda1, lambda2)

Lambda.list = list(Lambda1, Lambda2)

true.cluster = rbind(rep(1, n.subgroup), rep(1:2, each = n.subgroup/2))
true.beta <- beta
true.alpha <- alpha

# construct graph for the fused lasso
edge.matrix.full <- edge.generating(subgroup.name = 1:n.subgroup,
                                    vertex.size = 10)

main.simulation <- function (iter) {
  
  set.seed(iter)
  # generate data 
  data.simu <- simu.data.generating(n.subgroup, n.ober, n.T, n.X,
                                    Q.list, beta, alpha, error.sd, t.range = t.range)
  
  # oracle estimator --------------------------------------------------------
  t0 <- proc.time()
  fitted.oracle <- oracle.FLM(data = data.simu, Q.list = Q.list,
                              sp.basis = sp.basis)
  t.oracle <- proc.time() - t0 
  
  # individual group estimator ----------------------------------------------
  t0 <- proc.time()
  fitted.ind <- ind.FLM(data = data.simu, sp.basis)
  t.ind <- proc.time() - t0  
  
  # functional linear regression, k-means -----------------------------------
  t0 <- proc.time()
  fitted.kmeans <- kmeans.fusionFLM(data = data.simu, sp.basis)
  t.kmeans <- proc.time() - t0
  
  # functional linear regression, full graph --------------------------------
  t0 <- proc.time()
  fitted.full.kmeans <- fusionFLM(data = data.simu, sp.basis, edge.matrix.full,
                                  Lambda.list, initial.type = 'kmeans', 
                                  objective.path = FALSE, save.plot = FALSE)
  t.full.kmeans <- proc.time() - t0
  
  t0 <- proc.time()
  fitted.full.ind <- fusionFLM(data = data.simu, sp.basis, edge.matrix.full,
                               Lambda.list, initial.type = 'individual', 
                               objective.path = FALSE, save.plot = FALSE)
  t.full.ind <- proc.time() - t0
  
  t0 <- proc.time()
  fitted.full.lasso <- fusionFLM(data = data.simu, sp.basis, edge.matrix.full,
                                 Lambda.list, initial.type = 'lasso', 
                                 objective.path = FALSE, save.plot = FALSE)
  t.full.lasso <- proc.time() - t0
  
  
  # evaluation -------------------------------------------------------------
  result.ind <- eval.FLM(fitted.ind, true.cluster, true.beta, true.alpha, sp.basis)
  result.oracle <- eval.FLM(fitted.oracle, true.cluster, true.beta, true.alpha, sp.basis)
  result.kmeans <- eval.FLM(fitted.kmeans, true.cluster, true.beta, true.alpha, sp.basis)
  result.full.kmeans <- eval.FLM(fitted.full.kmeans, true.cluster, true.beta, true.alpha, sp.basis)
  result.full.ind <- eval.FLM(fitted.full.ind, true.cluster, true.beta, true.alpha, sp.basis)
  result.full.lasso <- eval.FLM(fitted.full.lasso, true.cluster, true.beta, true.alpha, sp.basis)
  
  time.all <- c(t.oracle[3], t.kmeans[3], t.ind[3], 
                t.full.kmeans[3], t.full.ind[3], t.full.lasso[3])
  
  cluster.all <- rbind(result.oracle$cluster.results,
                       result.ind$cluster.results,
                       result.kmeans$cluster.results,
                       result.full.kmeans$cluster.results,
                       result.full.ind$cluster.results,
                       result.full.lasso$cluster.results)
  
  rmse.beta.all <- rbind(result.oracle$rmse.beta,
                         result.ind$rmse.beta,
                         result.kmeans$rmse.beta,
                         result.full.kmeans$rmse.beta,
                         result.full.ind$rmse.beta,
                         result.full.lasso$rmse.beta)
  
  rmise.alpha.all <- rbind(result.oracle$rmise.alpha,
                           result.ind$rmise.alpha,
                           result.kmeans$rmise.alpha,
                           result.full.kmeans$rmise.alpha,
                           result.full.ind$rmise.alpha,
                           result.full.lasso$rmise.alpha)
  
  list(time.all, cluster.all, rmse.beta.all, rmise.alpha.all)
}
r.test = main.simulation(3)
print(r.test)
# save(file = 'result/test.RData', r.test)
# result <- mclapply(1:200, main.simulation, mc.cores = 16)
# save(file = paste0('result/cluster2_n', n.ober, 'sigma', error.sd, '.RData'), result)