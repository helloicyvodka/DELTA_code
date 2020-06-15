
mut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/mut-tree-tips",full.names = T)

wt.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/wt-tree-tips",full.names = T)




cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018_same.tsv"



#
outfile.parent.folder.mt_mt <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/same/mut-vs-mut/" 
outfile.parent.folder.mt_wt <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/same/mut-vs-wt/" 
outfile.parent.folder.wt_wt <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/same/wt-vs-wt/" 


parallel::mclapply(mut.genes, function(x) {
  treeS <- mut.files.list[grep(x, mut.files.list)]
  treeT <- wt.files.list[grep(x, wt.files.list)]
  outfile <- paste0(outfile.parent.folder.mt_mt, x, "_.almg")
  
  DELTA(
    treeS = treeS,
    treeT = treeS,
    cost = cost,
    outfile = outfile,
    method = "g",
    test = 10
  )
},
mc.cores = 20,
mc.silent = T)





parallel::mclapply(mut.genes, function(x) {
  treeS <- mut.files.list[grep(x, mut.files.list)]
  treeT <- wt.files.list[grep(x, wt.files.list)]
  outfile <- paste0(outfile.parent.folder.mt_wt, x, "_.almg")
  
  DELTA(
    treeS = treeS,
    treeT = treeT,
    cost = cost,
    outfile = outfile,
    method = "g",
    test = 10
  )
},
mc.cores = 20,
mc.silent = T)






parallel::mclapply(mut.genes, function(x) {
  treeS <- mut.files.list[grep(x, mut.files.list)]
  treeT <- wt.files.list[grep(x, wt.files.list)]
  outfile <- paste0(outfile.parent.folder.wt_wt, x, "_.almg")
  
  DELTA(
    treeS = treeT,
    treeT = treeT,
    cost = cost,
    outfile = outfile,
    method = "g",
    test = 10
  )
},
mc.cores = 20,
mc.silent = T)







# calculate mutated vs cut wt

mclapply(mut.files.list.genes,function(x){
  
  
  treeS <- mut.files.list[grep(x,mut.files.list)]
  treeT <- wt.cut.files.list[grep(x,wt.cut.files.list)]
  outfile <- paste0(outfile.parent.folder,x,"_.almg")

  DELTA(treeS=treeS,
        treeT=treeT,
        cost=cost,
        outfile=outfile,
        method = "g",
        test = 0
  )
  
  
},
mc.cores = 20,
mc.silent = T)



outfile.parent.folder <- "./data3/same/mut-vs-wt/"




result.list.files <- list.files( outfile.parent.folder,full.names = T)



result.list.files.genes <- unlist(sapply(list.files( outfile.parent.folder ,full.names = F),
                                         function(x){strsplit(x,split = "_")[[1]][1]}))


score.list <- parallel::mclapply(result.list.files,ggvita::readal.almg,mc.cores = 20)



names(score.list) <- result.list.files.genes



score.list.tip.size <- unlist(parallel::mclapply(wt.cut.files.list ,
                                                 function(x)nrow(as.data.frame(data.table::fread(x))),mc.cores = 20))



names(score.list.tip.size) <- sapply(wt.cut.files.list.genes,function(x)strsplit(x,split = "_")[[1]][1])




score.list.Score <- unlist( sapply(names(score.list),function(x){
  
  as.numeric(score.list[[x]]$Score)
  
}) )


score.df <- data.frame(gene=names(score.list),
                       DELTA.score=score.list.Score,
                       stringsAsFactors = F)


score.df$tip.size <- sapply(score.df$gene,function(x)score.list.tip.size[[x]])

score.df$score.divide.tip <- score.df$DELTA.score / score.df$tip.size

row.names(score.df) <- 1:nrow(score.df)





ID.df <- as.data.frame(data.table::fread("~/2017-2018/tree/tree_fate_change/code/id_transform/worm_proteinID_genename2.txt"))

ID.df$WikiGene_name <- sapply(ID.df$WikiGene_name,toupper)


score.df$gene_name <- sapply(score.df$gene ,function(x){
  
  if( x %in% ID.df$WikiGene_name){
        y <- ID.df[ID.df$WikiGene_name==x,"Protein_stable_ID"]
      }else{
        y <- x
      }
  return(y)
})


score.df.2 <- ID.df[ID.df$Protein_stable_ID %in% score.df$gene | ID.df$WikiGene_name %in% score.df$gene,]

score.df.2$DELTA.score <- sapply( 1:nrow(score.df.2), function(x){
  
  r <- score.df.2[x,]
  
  if(r[2] %in% score.df$gene ){ 
    y <- score.df[score.df$gene==as.character(r[2]),2]
  }
  else{
    y <- score.df[score.df$gene==as.character(r[1]),2]
  }
  
  return(y)

  })



score.df.2$tip.size <- sapply( 1:nrow(score.df.2), function(x){
  
  r <- score.df.2[x,]
  
  if(r[2] %in% score.df$gene ){ 
    y <- score.df[score.df$gene==as.character(r[2]),3]
  }
  else{
    y <- score.df[score.df$gene==as.character(r[1]),3]
  }
  
  return(y)
  
})


score.df.2$score.divide.tip <- sapply( 1:nrow(score.df.2), function(x){
  
  r <- score.df.2[x,]
  
  if(r[2] %in% score.df$gene ){ 
    y <- score.df[score.df$gene==as.character(r[2]),4]
  }
  else{
    y <- score.df[score.df$gene==as.character(r[1]),4]
  }
  
  return(y)
  
})




score.df.2$Protein_stable_ID.2 <- sapply( score.df.2$Protein_stable_ID, function(x){
  
   y <- strsplit(x,split = "\\.")[[1]]
   
  if(length(y)>2){
    y <- y[c(1,2)]
    y <- paste(y, collapse  = ".")
  }else{
      y <- x
  }
  
  y
  
})


toPlot.worm.2 <- toPlot.worm[toPlot.worm$orf %in% score.df.2$Protein_stable_ID.2,]

toPlot.worm.2$DELTA.score <- sapply(toPlot.worm.2$orf,function(x){
  
  unique(score.df.2[score.df.2$Protein_stable_ID.2==x,"DELTA.score"])
  
})

toPlot.worm.2$tip.size <- sapply(toPlot.worm.2$orf,function(x){
  
  unique(score.df.2[score.df.2$Protein_stable_ID.2==x,"tip.size"])
  
})


toPlot.worm.2$score.divide.tip <- sapply(toPlot.worm.2$orf,function(x){
  
  unique(score.df.2[score.df.2$Protein_stable_ID.2==x,"score.divide.tip"])
  
})








