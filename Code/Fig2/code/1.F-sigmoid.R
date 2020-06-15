F.sigmoid <- function(x, a=100){
  return (2 / (1 + exp(-a*x)) -1)
}