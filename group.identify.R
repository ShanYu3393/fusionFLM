group.identify <- function(eta, coefficient, group.length, edge.matrix) {
  
  start <- 1
  end <- 0
  cluster.num <- c()
  coeff.hat <- c()
  cluster.grp <- c()
  n.X <- length(group.length) - 1
  
  if(is.matrix(edge.matrix)) {
    tmp <- list()
    for(iter in 1:(1+n.X)) tmp[[iter]] <- edge.matrix
    edge.matrix <- tmp   
  }

  for (iter in 1:length(group.length)) {
    target <- apply(edge.matrix[[iter]], 1,
                    function(x) which(x == 1))
    source <- apply(edge.matrix[[iter]], 1,
                    function(x) which(x == -1))
    n.edge.matrix <- nrow(edge.matrix[[iter]])
    end <- n.edge.matrix * group.length[iter] + end
    eta.iter <- matrix(eta[start:end], ncol = n.edge.matrix)
    eta.norm <- apply(eta.iter, 2, function(x) sum(x^2))
    graph <- matrix(c(cbind(target, source)[which(eta.norm == 0), ]), ncol = 2)
    g <- graph.data.frame(graph, directed = F, vertices = 1:ncol(edge.matrix[[iter]]))
    Groups <- igraph::groups(components(g))
    CLUSTER <- rep(0, ncol(edge.matrix[[iter]]))
    for (igroup in 1:length(Groups)) {
      CLUSTER[as.numeric(Groups[[igroup]])] <- igroup
    }
    # calculate mean by group
    theta.matrix <- matrix(coefficient[[iter]], nrow = group.length[iter])
    coeff.mean <- aggregate(t(theta.matrix), list(CLUSTER), mean)[, -1]
    coeff.hat <- c(coeff.hat, t(as.matrix(coeff.mean)[CLUSTER, ]))
    cluster.num <- c(cluster.num, max(CLUSTER))
    cluster.grp <- rbind(cluster.grp, CLUSTER)
    cat(CLUSTER, "\n")
    start <- end + 1
  }
  list(cluster.num = cluster.num, coeff.hat = coeff.hat, cluster.grp = cluster.grp)
}
