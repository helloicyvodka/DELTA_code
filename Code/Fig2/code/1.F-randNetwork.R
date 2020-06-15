F.randNetwork <- function(N,K,digits = NULL){
  gene.interactions = rnorm(N*K)
  gene.network = c(gene.interactions, rep(0, N*(N-K)))
  if (is.null(digits))
    return (matrix(sample(gene.network), ncol = N, nrow = N))
  else
    return (round((matrix(sample(gene.network), ncol = N, nrow = N)), digits))
}

