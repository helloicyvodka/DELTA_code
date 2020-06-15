F.geneOF <- function(gene.level){
  if (any(gene.level > 1) || any(gene.level < -1)){
    print("warning! some gene.level is outside the value limit of -1<= gene.level <= 1")
  }
  return (as.numeric(gene.level > 0))
}
