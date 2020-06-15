F.treeCellType <- function(tree1){
  gene.OF <- lapply(tree1$leaves, function(node) paste0(node$Gene.OF,collapse = ""))
  return(length(unique(unlist(gene.OF))))
}
