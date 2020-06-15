library(plyr)
library(dplyr)
library(data.tree)
library(parallel)
library(rlist)




F.main <- function(ID,folder=".",
                  a,N,K,
                  Dmax, tmax, 
                  bias,RdropRate=0,
                  leafTypeMin,leafCountMin,tryMax,
                  digits=3,score="common",silent=F){
  
  
  
  Info  <-  paste0("ID", ID,
                      "a", a,
                      "N", N,
                      "K", K,
                      "t", tmax,
                      "D", Dmax,
                      "b", bias,
                      "R", RdropRate)
  
  parent.folder <- paste0(folder, "/",Info)
  
  dir.create(parent.folder)  

  outfilePrefix <- paste0(parent.folder,"/",Info)
  
  tr<- F.generateRandTreeWithFilter(a = a, N = N, K = K, Dmax = Dmax, 
                                     tmax = tmax, bias=bias, digits = digits,
                                     leafTypeMin = leafTypeMin,
                                     leafCountMin = leafCountMin, tryMax = tryMax)
  if(is.null(tr)){stop("No tree!")}
  
  cell.tree <- tr$cell.tree
  
  cell.tree.dfAll <- ToDataFrameTree(cell.tree,"pathString","level","name",'Lineage', 'DAge','tAge',
                                     Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                                     Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                                     isLeaf = function(node) node$isLeaf)


  
  cell.ReLin <- F.Simulate.Dropped.Leaf(cell.tree, RdropRate=RdropRate)
  
  outfile.ReLin.Leaf  <- paste0(outfilePrefix, ".ReLin.Leaf")
  
  outfile.ReLin.All  <- paste0(outfilePrefix, ".ReLin.All")
  
  write.table(cell.ReLin$List.ReLin$All,outfile.ReLin.All ,sep="\t",quote = F,col.names = T,row.names = F)
  
  write.table(cell.ReLin$List.ReLin$Leaf,outfile.ReLin.Leaf ,sep="\t",quote = F,col.names = T,row.names = F)
  

  cell.tree.Leaf <- cell.ReLin$File.alm
  
  cell.cost.file <- F.create.costfile(cell.tree.Leaf,score = score,silent = silent)
  
  
  outfile.network  <- paste0(outfilePrefix, ".cell.tree.network")
  write.table(tr$gene.network, outfile.network, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  outfile.cell.tree.all <- paste0(outfilePrefix, ".cell.tree.all")
  write.table(cell.tree.dfAll, outfile.cell.tree.all, sep="\t", row.names = FALSE)
  
  
  outfile.cell.tree.leaves <- paste0(outfilePrefix, ".alm")
  write.alm(cell.tree.Leaf,outfile.cell.tree.leaves)
  
  
  outfile.cell.cost.file <- paste0(outfilePrefix, ".tsv")
  write.table(cell.cost.file, outfile.cell.cost.file ,sep="\t",quote = F,col.names = F,row.names = F)
  
  
  
  fileInfo <- paste0(parent.folder,"/",Info,".info")
  
  fileInfoWrite <- file(fileInfo,"w")
 
  

  
  writeLines(text = c(paste0("ID=",ID),
                      paste0("a=",a),
                      paste0("N=",N),
                      paste0("K=",K),
                      paste0("tmax=",tmax),
                      paste0("Dmax=",Dmax),
                      paste0("bias=",bias),
                      paste0("RdropRate=",RdropRate),
                      paste0("Score type=",score),
                      paste0("LeafCount=",tr$leafCount),
                      paste0("LeafTypeCount=",tr$leafTypeCount),
                      "\n",
                      "gene.int :",
                      paste(tr$gene.init,collapse = " "),
                      "\n",
                      "gene.network :",
                      unlist(apply(tr$gene.network,1,function(x)paste(as.vector(x),collapse ="\t")))),
             con = fileInfoWrite)
             
  close(fileInfoWrite)
  
  if(!silent) print(paste0("write the result to file ", normalizePath(parent.folder)))

  
}