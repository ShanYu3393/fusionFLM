jaccard <- function(group1, group2) {
  n <- length(group1)
  result <- c()
  if (length(group1) == length(group2)) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (group1[i] == group1[j] & group2[i] == group2[j]) {
          result <- c(result, 1)
        } else if (group1[i] != group1[j] & group2[i] != group2[j]) {
          result <- c(result, -1)
        } else {
          result <- c(result, 0)
        }
      }
    }
  }
  sum(result > 0) / sum(result >= 0)
}

# jaccard(group1=c(1,1,1,2,2),group2=c(2,2,2,3,3))
