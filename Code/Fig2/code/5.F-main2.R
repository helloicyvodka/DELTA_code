library(plyr)
library(dplyr)
library(data.tree)
library(parallel)
library(rlist)




F.main.2 <- function(ID,folder=".",
                     a,N,K,
                     Dmax, tmax, 
                     bias,RdropRate=0,
                     leafTypeMin,leafCountMin,tryMax,
                     digits=3,score="common",m,u,silent=F){
  

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
  
  tr <- F.generateRandTreeWithFilter.2(a = a, N = N, K = K, Dmax = Dmax, 
                                    tmax = tmax, bias=bias, digits = digits,
                                    leafTypeMin = leafTypeMin,
                                    leafCountMin = leafCountMin, tryMax = tryMax)
  

  
  if(is.null(tr)){stop("No tree!")}
  
  cell.tree.1 <- tr$cell.tree$cell.tree.1
  cell.tree.2 <- tr$cell.tree$cell.tree.2
  
  cell.tree.dfAll.1 <- ToDataFrameTree(cell.tree.1,"pathString","level","name",'Lineage', 'DAge','tAge',
                                     Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                                     Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                                     isLeaf = function(node) node$isLeaf)
  
  cell.tree.dfAll.2 <- ToDataFrameTree(cell.tree.2,"pathString","level","name",'Lineage', 'DAge','tAge',
                                     Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                                     Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                                     isLeaf = function(node) node$isLeaf)
  
  
  cell.ReLin.1 <- F.Simulate.Dropped.Leaf(cell.tree.1, RdropRate=RdropRate)
  cell.ReLin.2 <- F.Simulate.Dropped.Leaf(cell.tree.2, RdropRate=RdropRate)
  
  outfile.ReLin.Leaf.1  <- paste0(outfilePrefix, ".1.ReLin.Leaf")
  outfile.ReLin.Leaf.2  <- paste0(outfilePrefix, ".2.ReLin.Leaf")
  
  outfile.ReLin.All.1  <- paste0(outfilePrefix, ".1.ReLin.All")
  outfile.ReLin.All.2  <- paste0(outfilePrefix, ".2.ReLin.All")
  
  write.table(cell.ReLin.1$List.ReLin$"All",outfile.ReLin.All.1 ,sep="\t",quote = F,col.names = T,row.names = F)
  write.table(cell.ReLin.2$List.ReLin$"All",outfile.ReLin.All.2 ,sep="\t",quote = F,col.names = T,row.names = F)
  
  write.table(cell.ReLin.1$List.ReLin$"Leaf",outfile.ReLin.Leaf.1 ,sep="\t",quote = F,col.names = T,row.names = F)
  write.table(cell.ReLin.2$List.ReLin$"Leaf",outfile.ReLin.Leaf.2 ,sep="\t",quote = F,col.names = T,row.names = F)
  
  
  cell.tree.Leaf.1 <- cell.ReLin.1$File.alm
  cell.tree.Leaf.2 <- cell.ReLin.2$File.alm
  cell.tree.Leaf   <- rbind(cell.tree.Leaf.1,cell.tree.Leaf.2)

  
  outfile.network  <- paste0(outfilePrefix, ".cell.tree.network")
  write.table(tr$gene.network, outfile.network, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  outfile.cell.tree.all.1 <- paste0(outfilePrefix, ".1.cell.tree.all")
  write.table(cell.tree.dfAll.1, outfile.cell.tree.all.1, sep="\t", row.names = FALSE)
  outfile.cell.tree.all.2 <- paste0(outfilePrefix, ".2.cell.tree.all")
  write.table(cell.tree.dfAll.2, outfile.cell.tree.all.2, sep="\t", row.names = FALSE)
  
  outfile.cell.tree.leaves.1 <- paste0(outfilePrefix, ".1.alm")
  write.alm(cell.tree.Leaf.1,outfile.cell.tree.leaves.1)
  
  outfile.cell.tree.leaves.2 <- paste0(outfilePrefix, ".2.alm")
  write.alm(cell.tree.Leaf.2,outfile.cell.tree.leaves.2)
  
  # m: match loci score
  # u: unmatch loci score
  
  cell.cost.file <- F.create.costfile(cell.tree.Leaf,score = score,silent = silent,m=m,u=u)
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
                      paste0("LeafCountMin=",tr$leafCount),
                      paste0("LeafTypeCountMin=",tr$leafTypeCount),
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
