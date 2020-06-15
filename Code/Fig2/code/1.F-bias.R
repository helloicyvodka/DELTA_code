
F.bias <- function(N,bias){
  sapply(runif(N-2,0,1),function(x)ifelse(x>=bias,1,-1))
}
