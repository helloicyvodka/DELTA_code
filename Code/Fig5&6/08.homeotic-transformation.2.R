
mut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/mutated_for_DELTA/mutated_for_DELTA",full.names = T)

mut.files.list.genes <- unlist( sapply(list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/mutated_for_DELTA/mutated_for_DELTA",full.names = F),function(x){strsplit(x,split = "_")[[1]][1]}) )

wt.cut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/wt_cut_for_DELTA",full.names = T)

wt.cut.files.list.genes <- unlist( sapply( list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/wt_cut_for_DELTA",full.names = F),function(x){strsplit(x,split = "_")[[1]][1]}) )

wt.file <- "~/2017-2018/tree/tree_fate_change/data/data_from_digital/WT_merge/wt_for_DELTA.alm"

cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"

cost_mt10_unmt1 <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost_mt10_unmt1.tsv"

cost_mt10_unmt0<- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost_mt10_unmt0.tsv"



hom.trans.outfile.parent.folder <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/results/DELTA_results/20180528_homeotic_transformation_score" 


create.cost <- function(mt.s, unmt.s,file){
  
  
  tr.Class <- c(0:7)
  
  cost.mt <-as.data.frame(expand.grid(tr.Class,tr.Class))
  
  
  Dec2Bin <- function(x) {
    i <- 0
    string <- numeric(32)
    if(x==0){return(x)}
    while(x > 0) {
      string[32 - i] <- x %% 2
      x <- x %/% 2
      i <- i + 1
    }
    first <- match(1, string)
    paste(as.character(string[first:32]),collapse="")
  }
  
  
  
  to.bin.vec <- function(x){
    
    x <- Dec2Bin(x)
    
    if(nchar(x) == 1){
      x <- as.vector(as.numeric(x))
    }else{
      x <- as.numeric(unlist(strsplit(x, split = "")))
    }
    
    
    if( length(x) <3 ){
      x <- c(rep(0,3-length(x)),x)
    }
    
    return(x)
  }
  
  
  
  
  
  
  cost.mt$score <- sapply(1:nrow(cost.mt),function(x){
    
    r <- cost.mt[x,]
    
    r.1 <- to.bin.vec(r[1])
    r.2 <- to.bin.vec(r[2])
    
    r.r <- abs(r.1-r.2)
    
    y <- length(r.r[r.r==0])*mt.s+length(r.r[r.r==1])*(-unmt.s)
    
    return(y)
    
  })
  
  write.table(as.matrix(cost.mt),file = file ,quote = F,row.names = F,col.names = F)
  
  
}

cost  <- "./data/data_from_digital/cost_mt10_unmt1.tsv"

create.cost(10,1,cost)

cal.S_trans <- function(x,prune=1,test=0,max_target=1000,cost,mc.cores=5,...){
  
  x <- toupper(x)
  
  tree.mt <- mut.files.list[grep(x,mut.files.list)]
  tree.wt <- wt.cut.files.list[grep(x,wt.cut.files.list)]
  outfile.mt_wt <- paste0(hom.trans.outfile.parent.folder,"/mt_wt/",x,"_.alml")
  outfile.mt_mt <- paste0(hom.trans.outfile.parent.folder,"/mt_mt/",x,"_.alml")
  outfile.wt_wt <- paste0(hom.trans.outfile.parent.folder,"/wt_wt/",x,"_.alml")
  
  l <- list(list(tree.wt,tree.wt,outfile.wt_wt),list(tree.mt,tree.wt,outfile.mt_wt),list(tree.mt,tree.mt,outfile.mt_mt))
  
  mclapply(l,function(x){
    
    DELTA(treeS=x[[1]],
          treeT=x[[2]],
          cost=cost,
          outfile=x[[3]],
          method = "l",
          prune = prune,
          test = test,
          max_target = max_target,
          ...
    )
    
  },mc.cores = mc.cores)

  hom.trans.S_trans <- list()
  hom.trans.S_trans[[x]] <- mclapply(list(outfile.wt_wt,
                                          outfile.mt_wt,
                                          outfile.mt_mt),
                                     readal.alml2,
                                     mc.cores = 5)
  
  names(hom.trans.S_trans[[x]]) <- c("wt_wt","mt_wt","mt_mt")
  
  
  
  hom.trans.S_trans[[x]]$digital.pairs <- 
    hom.trans.dt %>% 
    filter(toupper(Gene)==x) %>% 
    list.parse()
  
  
  
  hom.trans.S_trans[[x]]$S.wt_wt <-
    hom.trans.S_trans[[x]]$digital.pairs %>% 
    lapply(., function(i){
      hom.trans.S_trans[[x]]$wt_wt %>% 
        list.filter(
          RootS==i$from_cell_LN & RootT==i$to_cell_LN) %>% 
        list.select(Score,RootT,RootS) 
    }) 
  
  
  
  
  hom.trans.S_trans[[x]]$S.mt_wt <-
    hom.trans.S_trans[[x]]$digital.pairs %>% 
    lapply(., function(i){
      hom.trans.S_trans[[x]]$mt_wt %>% 
        list.filter(
          RootS==i$from_cell_LN & RootT==i$to_cell_LN) %>% 
        list.select(Score,RootT,RootS) 
    })
  
  
  hom.trans.S_trans[[x]]$S.mt_mt <-
    hom.trans.S_trans[[x]]$digital.pairs %>% 
    lapply(., function(i){
      hom.trans.S_trans[[x]]$mt_mt %>% 
        list.filter(
          RootS==i$from_cell_LN & RootT==i$to_cell_LN) %>% 
        list.select(Score,RootT,RootS) 
    })
  
  
  return(hom.trans.S_trans)
  
}



wt.cut.files.list.genes <- unlist( sapply( list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/wt_cut_for_DELTA",full.names = F),function(x){strsplit(x,split = "_")[[1]][1]}) )

hom.trans.gene.avalibale <- intersect(toupper(as.character(unique(hom.trans.dt$Gene))),toupper(wt.cut.files.list.genes)) %>% setdiff(.,c("CIR-1","D1054.14","F53B7.3","HIM-10","R144.2","SRC-1"))

hom.trans.S_trans.1000 <-
  mclapply(hom.trans.gene.avalibale,
           function(x)
             cal.S_trans(
               x,
               prune=5,
               max_target = 1000,
               mc.cores = 3,
               cost = cost
             ),
           mc.cores = 12)


names(hom.trans.S_trans.1000) <- as.character(hom.trans.gene.avalibale)

hom.trans.dt.ava <- hom.trans.dt %>% filter(toupper(Gene) %in% hom.trans.gene.avalibale)

hom.trans.dt.ava$S.wt_wt <- sapply(1:nrow(hom.trans.dt.ava),function(x){
  
 r <- hom.trans.dt.ava[x,]
  
 S.wt_wt <-  hom.trans.S_trans.1000[[toupper(as.character(r[,1]))]][[toupper(as.character(r[,1]))]]$S.wt_wt

 s <- sapply(S.wt_wt,function(z){
   
   z %>%
     list.filter(     RootS == as.character(r[,4]) & 
                      RootT == as.character(r[,5]) ) %>%
     list.select(Score) %>%
     unlist() %>%
     as.numeric() %>%
     unique()
   
  }) %>% 
   unlist() %>% 
   unique()
 
   return(ifelse(length(s)!=0,s,NA))
 
})

hom.trans.dt.ava$S.mt_wt <- sapply(1:nrow(hom.trans.dt.ava),function(x){
  
  r <- hom.trans.dt.ava[x,]

  S.mt_wt <-  hom.trans.S_trans.1000[[toupper(as.character(r[,1]))]][[toupper(as.character(r[,1]))]]$S.mt_wt

  s <- sapply(S.mt_wt,function(z){
    
    z %>%
      list.filter(     RootS == as.character(r[,4]) & 
                       RootT == as.character(r[,5]) ) %>%
      list.select(Score) %>%
      unlist() %>%
      as.numeric() %>%
      mean()
    
  }) %>% 
    unlist() %>% 
    unique()
  
  return(ifelse(length(s)!=0,s,NA))
  
})

hom.trans.dt.ava$S.mt_mt <- sapply(1:nrow(hom.trans.dt.ava),function(x){
  
  r <- hom.trans.dt.ava[x,]
  
  S.mt_mt <-  hom.trans.S_trans.1000[[toupper(as.character(r[,1]))]][[toupper(as.character(r[,1]))]]$S.mt_mt
  
  s <- sapply(S.mt_mt,function(z){
    
    z %>%
      list.filter(     RootS == as.character(r[,4]) & 
                       RootT == as.character(r[,5]) ) %>%
      list.select(Score) %>%
      unlist() %>%
      as.numeric() %>%
      unique()
    
  }) %>% 
    unlist() %>% 
    unique()
  
  return(ifelse(length(s)!=0,s,NA))
  
})

hom.trans.dt.ava <- hom.trans.dt.ava %>% mutate(S.diff=S.mt_wt-(S.wt_wt+S.mt_mt)/2)

hom.trans.dt.ava <- hom.trans.dt.ava %>% mutate(S.ratio=S.mt_wt/((S.wt_wt+S.mt_mt)/2))

hom.trans.dt.ava %>% na.omit() %>% View()


sym.sisters <- read.csv("../transform ppp to 111/symmetry-sisters.csv",header=T,stringsAsFactors = F ,colClasses = "character")
sym.sisters <- data.frame(Gene =rep("Sym.sisters",nrow(sym.sisters)),sym.sisters)
sym.sisters.2 <- data.frame(sym.sisters[,c(1)],sym.sisters[,c(3)],sym.sisters[,c(2)],sym.sisters[,c(5)],sym.sisters[,c(4)])
colnames(sym.sisters.2) <- colnames(sym.sisters)
sym.sisters <- rbind(sym.sisters,sym.sisters.2)
colnames(sym.sisters) <- colnames(hom.trans.dt.ava)[1:5]
hom.and.sym.dt <- rbind(hom.trans.dt,sym.sisters)















