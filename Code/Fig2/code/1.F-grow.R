
F.grow <- function(gene.level, gene.network, a = 100,bias = bias){
  return (F.sigmoid(as.vector(gene.network %*% gene.level)*c(1,1,F.bias(length(gene.level),bias)), a))
}
