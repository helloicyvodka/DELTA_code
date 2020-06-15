F.generateRandTreeWithFilter  <- function(a, N, K, Dmax, tmax, bias, digits = 3, 
                                          leafTypeMin = 6, leafCountMin = 500, tryMax = 100){
  for (i in 1:tryMax){
    celltree <- F.generateRandTree(a = a, N = N, K=K, Dmax = Dmax, tmax = tmax, bias = bias, digits = digits)
    leafTypeCount <- F.treeCellType(celltree[["cell.tree"]])
    leafCount <- celltree[["cell.tree"]]$leafCount
    celltree$leafTypeCount <- leafTypeCount
    celltree$leafCount <- leafCount
    
    if (leafCount >= leafCountMin && leafTypeCount >= leafTypeMin){
      return(celltree)
    }
  }
  print(paste0("tried ", tryMax, " times. Still No good. return NULL"))
  return(NULL)
}




#N=12


F.generateRandTreeWithFilter.2  <- function(a, N, K, Dmax, tmax, bias, digits = 3, 
                                          leafTypeMin = 6, leafCountMin = 500, tryMax = 100){
  for (i in 1:tryMax){
    
    celltree        <- F.generateRandTree.2(a = a, N = N, K=K, Dmax = Dmax, tmax = tmax, bias = bias, digits = digits)
    celltree.1      <- celltree$cell.tree$cell.tree.1
    celltree.2      <- celltree$cell.tree$cell.tree.2

    leafTypeCount.1 <- F.treeCellType(celltree.1)
    leafTypeCount.2 <- F.treeCellType(celltree.2)
    leafTypeCount   <- min(leafTypeCount.1,leafTypeCount.2)
    
    leafCount.1     <- celltree.1$leafCount
    leafCount.2     <- celltree.2$leafCount
    leafCount       <- min(leafCount.1,leafCount.2)
    
    if (leafCount >= leafCountMin && leafTypeCount >= leafTypeMin){
      celltree$leafTypeCount <- leafTypeCount
      celltree$leafCount <- leafCount
      return(celltree)
    }
  }
  print(paste0("tried ", tryMax, " times. Still No good. return NULL"))
  return(NULL)
}







