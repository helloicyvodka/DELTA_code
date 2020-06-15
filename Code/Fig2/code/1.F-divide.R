F.divide <- function(gene.level){
  if (gene.level[1] > 0){
    L <-  gene.level
    R <- gene.level
    L[1] <- -1
    R[1] <- -1
    if (gene.level[2] >0) L[2] <- -1
    return (data.frame(L,R))
  }
  return (NULL)
}
