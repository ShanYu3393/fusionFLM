edge.generating <- function(subgroup.name, type = "full", knn = NULL,
                        dist = NULL, Group.est = 1, Group.true = '',
                        vertex.size = 6, cutoff.value = NULL) {

  if (type == "full") {
    
    delta.graph <- delta.matrix(length(subgroup.name))
    
    # generate delta matrix for full-connected graph
    g <- make_full_graph(length(subgroup.name))
    
  } else {
    if (is.null(dist)) stop("lack of distance matrix")

    if (type == "knn") {
      if (is.null(knn)) {
        stop("lack of knn")
      } else {
        KNN <- list()
        for (iter in 1:length(subgroup.name)) KNN[[iter]] <- order(dist[iter, ], decreasing = TRUE)[2:knn]

        graph.knn <- data.frame(source = rep(subgroup.name, each = knn - 1))
        target <- c()
        for (iter in 1:length(subgroup.name)) {
          target <- c(target, subgroup.name[KNN[[iter]]])
        } 
        graph.knn$target <- target

        A <- apply(graph.knn, 1, function(x) paste0(sort(x), collapse = "#"))
        B <- unique(A)
        graph.knn0 <- unlist(strsplit(B, "#"))
        graph.knn2 <- data.frame(
          source = graph.knn0[seq(1, length(graph.knn0) - 1, 2)],
          target = graph.knn0[seq(2, length(graph.knn0), 2)]
        )
        graph.knn2 <- graph.knn2[which(as.character(graph.knn2[, 1]) != as.character(graph.knn2[, 2])), ]

        vertices <- data.frame("name" = subgroup.name)
        g <- graph.data.frame(graph.knn2, directed = FALSE, vertices = vertices)
        vertices$group <- edge.betweenness.community(g)$membership
      }
      
      # generate delta matrix for knn graph
      delta.graph <- matrix(0, ncol = length(subgroup.name),
                            nrow = nrow(graph.knn2))
      for (iter in 1:nrow(graph.knn2)) {
        delta.graph[iter, which(graph.knn2[iter, 1] == subgroup.name)] <- 1
        delta.graph[iter, which(graph.knn2[iter, 2] == subgroup.name)] <- -1
      }
    }

    if (type == "mst") {
      A.tree <- ape::mst(exp(-dist))

      index <- unique(t(apply(which(A.tree != 0, arr.ind = T), 1, sort)))
      delta.graph <- c()
      for (iter in 1:nrow(index)) {
        row.iter <- rep(0, length(subgroup.name))
        row.iter[index[iter, 1]] <- 1
        row.iter[index[iter, 2]] <- -1
        delta.graph <- rbind(delta.graph, row.iter)
      }
      
      # generate delta matrix for mst tree
      g <- graph_from_edgelist(index, directed = FALSE)
  
    }
    
    if (type == "threshhold") {
      if (is.null(cutoff.value)) {
        stop("lack of cutoff value")
      } else {
        Neighbor <- list()
        for (iter in 1:length(subgroup.name)) {
          Neighbor[[iter]] <- which(dist[iter, ] > cutoff.value &
                                    dist[iter, ] != max(dist[iter, ]))
        } 
        
        graph.neighbor <- data.frame(source = rep(subgroup.name,
                                times = unlist(lapply(Neighbor, length))))
        target <- c()
        for (iter in 1:length(subgroup.name)) {
          target <- c(target, subgroup.name[Neighbor[[iter]]])
        } 
        graph.neighbor$target <- target
        
        A <- apply(graph.neighbor,
                   1, function(x) paste0(sort(x), collapse = "#"))
        B <- unique(A)
        graph.neighbor0 <- unlist(strsplit(B, "#"))
        graph.neighbor2 <- data.frame(
          source = graph.neighbor0[seq(1, length(graph.neighbor0) - 1, 2)],
          target = graph.neighbor0[seq(2, length(graph.neighbor0), 2)]
        )
        graph.graph.neighbor2 <- graph.neighbor2[which(as.character(graph.neighbor2[, 1]) != as.character(graph.neighbor2[, 2])), ]
        
        vertices <- data.frame("name" = subgroup.name)
        g <- graph.data.frame(graph.neighbor2, directed = FALSE, vertices = vertices)
      }
      
      # generate delta matrix for knn graph
      delta.graph <- matrix(0, ncol = length(subgroup.name),
                            nrow = nrow(graph.neighbor2))
      for (iter in 1:nrow(graph.neighbor2)) {
        delta.graph[iter, which(graph.neighbor2[iter, 1] == subgroup.name)] <- 1
        delta.graph[iter, which(graph.neighbor2[iter, 2] == subgroup.name)] <- -1
      }
    }
  }
  
  plot(g, vertex.label = Group.true,
       vertex.color = Group.est,
       vertex.size = vertex.size,
       vertex.label.cex = 0.8,
       edge.arrow.size = 0.8)

  return(delta.graph)
}
