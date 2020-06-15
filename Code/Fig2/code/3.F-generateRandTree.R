F.generateRandTree <- function(a, N, K, Dmax, tmax, bias, digits = 3){
  gene.init  <- F.randInit(N)
  gene.network <- F.randNetwork(N, K, digits = digits)
  cell.tree <- F.develop(gene.init, gene.network, tmax = tmax, Dmax = Dmax, a = a, bias = bias)
  return(list(gene.init = gene.init, gene.network = gene.network, cell.tree = cell.tree))
}



F.generateRandTree.2 <- function(a, N, K, Dmax, tmax, bias, digits = 3){
  gene.init  <- F.randInit(N)
  gene.network <- F.randNetwork(N, K, digits = digits)
  cell.tree.1 <- F.develop(gene.init, gene.network, tmax = tmax, Dmax = Dmax, a = a, bias = bias)
  cell.tree.2 <- F.develop(gene.init, gene.network, tmax = tmax, Dmax = Dmax, a = a, bias = bias)
  return(list(gene.init = gene.init, gene.network = gene.network, cell.tree = list(cell.tree.1=cell.tree.1,
                                                                                   cell.tree.2=cell.tree.2)))
}
