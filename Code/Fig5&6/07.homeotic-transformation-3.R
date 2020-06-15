library(rlist)
library(plyr)
library(dplyr)
library(parallel)

hom.trans.dt <- read.table("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/cell_to_cell_fateChange_data/melt_fateChange.txt",
                           header = T)
hom.trans.uni.dt <- hom.trans.dt[,c(2,3,4,5)] %>% list.parse() %>% unique() %>% list.stack()

hom.trans.uni.dt <- union(hom.trans.uni.dt[,c(3,4)] %>% list.parse(),hom.trans.uni.dt[,c(4,3)] %>% list.parse()) %>% unique() %>% list.stack()

hom.trans.uni.dt<- hom.trans.uni.dt %>% list.parse() %>% unique() %>% list.stack()

mut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/mutated_for_DELTA/mutated_for_DELTA",full.names = T)

mut.files.list.genes <- unlist( sapply(list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/mutated_for_DELTA/mutated_for_DELTA",full.names = F),function(x){strsplit(x,split = "_")[[1]][1]}) )

wt.cut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/wt_cut_for_DELTA",full.names = T)

wt.cut.files.list.genes <- unlist( sapply( list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/wt_cut_for_DELTA",full.names = F),function(x){strsplit(x,split = "_")[[1]][1]}) )

wt.file <- "~/2017-2018/tree/tree_fate_change/data/data_from_digital/WT_merge/wt_for_DELTA.alm"

cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"

hom.trans.outfile.parent.folder <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/results/DELTA_results/20180525_homeotic_transformation_rank/" 



mclapply(mut.files.list.genes,function(x){
  
  treeS <- mut.files.list[grep(x,mut.files.list)]
  treeT <- wt.cut.files.list[grep(x,wt.cut.files.list)]
  outfile <- paste0(hom.trans.outfile.parent.folder,x,"_.alml")
  
  DELTA(treeS=treeS,
        treeT=treeT,
        cost=cost,
        outfile=outfile,
        method = "l",
        prune = 1,
        test = 100,
        max_target = 1000
  )
  
  
},
mc.cores = 20,
mc.silent = T)




hom.trans.result.list.files <- list.files( hom.trans.outfile.parent.folder, full.names = T)

hom.trans.result.list.files.genes <- unlist(sapply(list.files( hom.trans.outfile.parent.folder,full.names = F),
                                         function(x){strsplit(x,split = "_")[[1]][1]}))




hom.trans.score.list <- list()

# in the next process, I delete some genes that can not DELTA out result.
hom.trans.fail.DELTA.gene <- c("CIR-1","D1054.14","F53B7.3","HIM-10","R144.2","SRC-1")   

for(i in 1:198){
  
  cat(i,"\n")
  hom.trans.score.list[[hom.trans.result.list.files.genes[i]]] <- readal.alml2(hom.trans.result.list.files[i])
  
}


hom.trans.dt.gene <- hom.trans.dt$Gene %>% unique()

hom.trans.rank.dt <- data.frame(Gene=as.character(hom.trans.dt.gene),Gene.upper=as.character(hom.trans.dt.gene) %>% toupper(),stringsAsFactors = F)

hom.trans.rank.dt$mean.rank <- sapply(hom.trans.rank.dt$Gene,function(i){
  
  dt <- hom.trans.dt %>% filter(Gene==i)
  
  l <- hom.trans.score.list[[toupper(i)]]
  
  rank.seq <- sapply(1:length(l),function(ii){
    
    ll <- l[[ii]]
    
    f.1 <- dt[dt$from_cell_LN==ll$RootS & dt$to_cell_LN==ll$RootT, ]
    f.2 <- dt[dt$from_cell_LN==ll$RootT & dt$to_cell_LN==ll$RootS, ]
    
    if(nrow(f.1)==1 | nrow(f.2)==1){
      
      y <- ii
      
    }else{ y <- NA }
    
    return(y)

  })
  
  yy <- mean(na.omit(rank.seq))
  
  return(yy)
  
  
  
})

#write.table(hom.trans.rank.dt,file="./data/homeotic_trans_rank/max_target_1000_pr_1.txt",quote = F,row.names = F,col.names = T)

hom.trans.rank.dt.omit <- hom.trans.rank.dt %>% na.omit()


hom.trans.rank.dt.omit$gene_name <- sapply( hom.trans.rank.dt.omit$Gene.upper,  function(x){
  
  if( x %in% ID.df$WikiGene_name){
    y <- ID.df[ID.df$WikiGene_name==x,"Protein_stable_ID"]
  }else{
    y <- x
  }

  return(y)
  
})


hom.trans.rank.dt.omit.2 <- ID.df[ID.df$Protein_stable_ID %in% hom.trans.rank.dt.omit$gene | ID.df$WikiGene_name %in% hom.trans.rank.dt.omit$gene,]

hom.trans.rank.dt.omit.2$DELTA.rank.mean <- sapply( 1:nrow(hom.trans.rank.dt.omit.2), function(x){
  
  r <- hom.trans.rank.dt.omit.2[x,]
  
  if(r[2] %in% hom.trans.rank.dt.omit$gene ){ 
    y <- hom.trans.rank.dt.omit[hom.trans.rank.dt.omit$gene==as.character(r[2]),3]
  }
  else{
    y <- hom.trans.rank.dt.omit[hom.trans.rank.dt.omit$gene==as.character(r[1]),3]
  }
  
  return(y)
  
})

hom.trans.toPlot <- merge(hom.trans.rank.dt.omit.2,toPlot.worm.2,by.x="Protein_stable_ID",by.y="orf")



ggplot(hom.trans.toPlot,aes(x=dN,y=DELTA.rank.mean))+geom_point()+geom_smooth(method = "loess")




data.for.plot <- hom.trans.toPlot

panel.cor<- function(x, y, digits = 2, cex.cor, ...)
{
  df <- data.frame(x=x,y=y)
  df <- na.omit(df)
  x <- df$x
  y <- df$y
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  r <- abs(cor(x, y,method = "spearman"))
  txt.r <- format(c(r, 0.123456789), digits = digits)[1]
  txt.r <- paste0("R= ", txt.r)
  
  p <- cor.test(x, y,method = "spearman")$p.value
  txt.p <- format(c(p, 0.123456789), digits = digits)[1]
  if(p<0.05) txt.p <- paste0("p.value ","\n","< 0.05")
  else txt.p <- paste0("p.value= ","\n", txt.p)
  
  if(missing(cex.cor)) cex.cor.r <- 0.8/strwidth(txt.r);cex.cor.p <- 0.8/strwidth(txt.p)
  
  text(0.5, 0.8, txt.r, cex = 1)
  text(0.5,0.2,txt.p, cex = 1)
}





panel.abline <-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = 2, col = col, bg = bg, cex = 0.1)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y~x), 
           col = col.smooth, ...)
}



panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

pairs(data.for.plot[,c(3,4,8:13)],
      lower.panel = panel.smooth, 
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      cex.labels = 0.8,
      main="Mutated Tree")






