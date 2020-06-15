F.randInit <- function(N){
  randInit <- c(-1,sample(c(1,-1), size = N-1, replace = TRUE))
  return (randInit)
}