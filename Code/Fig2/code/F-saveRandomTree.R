F.saveRandomTree <- function(cell.tree, network = NULL, score = "common", 
                             class = "Gene.OF", outfilePrefix = "",silent = FALSE) {
  
  # transform random tree to a data.frame
  cell.tree.dfAll <- ToDataFrameTree(cell.tree,"pathString","level","name",'Lineage', 'DAge','tAge',
                                     Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                                     Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                                     isLeaf = function(node) node$isLeaf)
  
  outfile.cell.tree.all <- paste0(outfilePrefix, "cell.tree.all")
  
  # save the data.frame as a file
  write.table(cell.tree.dfAll, outfile.cell.tree.all, sep="\t", row.names = FALSE)
  
  if (!silent) print(paste0("write the all details for the tree to file ", outfile.cell.tree.all))
  
  cell.tree.dfLeaf <- ToDataFrameTable(cell.tree,"pathString","level","name",'Lineage', 'DAge','tAge',
                                       Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))
  
  cell.tree.Leaf <- cell.tree.dfLeaf[,c("Lineage","name","Gene.OF")]
  
  outfile.cell.tree.leaves <- paste0(outfilePrefix, "cell.tree.leaves")
  
  write.table(cell.tree.Leaf, outfile.cell.tree.leaves,sep = "\t", row.names = FALSE, quote = FALSE)
  
  #remove the last newline symbol in outfile.cell.tree.leaves
  
  txt <- readChar(outfile.cell.tree.leaves,file.info(outfile.cell.tree.leaves)$size,useByte = TRUE)
  
  txt <- gsub("\r","",txt)
  
  if (substr(txt,nchar(txt),nchar(txt)) == "\n")
    cat(substr(txt,1,nchar(txt)-1), file = outfile.cell.tree.leaves)
  
  if (!silent) print(paste0("write the leaves file for HSA to ", outfile.cell.tree.leaves))
  
  outfile.cell.tree.score <- paste0(outfilePrefix, "cell.tree.score")
  
  gene.OF <- cell.tree.dfLeaf$Gene.OF
  
  if(length(unique(gene.OF)) > 2){
    gene.OF.pairs <- combn(unique(gene.OF),2)
  }else{
    gene.OF.pairs <- NULL
  }
  
  
  dfscore <- data.frame(S=c(unique(gene.OF),gene.OF.pairs[1,]), T=c(unique(gene.OF),gene.OF.pairs[2,]))
  
  
  if (score == "common"){
    # Hamming distance 
    dfscore$Score <- apply(dfscore,1,function(x) sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],""))))
  }
  else if (score == "same"){
    dfscore$Score <- (dfscore$S == dfscore$T)*2
  }
  else if (score == "norm"){
    dfscore$Score <- apply(dfscore,1,function(x) sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],""))))
    dfscore$Score <- round((dfscore$Score - mean(dfscore$Score))/sd(dfscore$Score)) + 1
  } 
  else if (!silent) print("invalid value for score, score can be same, common or norm")
  
  write.table(dfscore, outfile.cell.tree.score, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  if (!silent) print(paste0("write the score file to ", outfile.cell.tree.score))
  
  if (!is.null(network)) {
    outfile.network  <- paste0(outfilePrefix, "cell.tree.network")
    write.table(network, outfile.network, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste0("write the interaction network to file ", outfile.network))
  }else if (!silent) print("You need to provide Network in order to save it!")
  
  
  
}
