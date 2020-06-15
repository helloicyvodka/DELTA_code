

# most usual used 

setwd("~/2017-2018/tree/tree_fate_change")

library(plyr)
library(dplyr)
library(ggplot2)
library(ggvita)
library(parallel)

mut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/mut-tree-tips",full.names = T)

wt.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/wt-tree-tips",full.names = T)

cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"



# data are in  "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/" (tab to find)
hom.trans.outfile.parent.folder <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/common-alml/"
prune=1
Cal=F


result.list.files.genes  <- unlist(sapply(list.files("./data3/same-almg/mut-vs-mut/",full.names = F),
                                                                           function(x){strsplit(x,split = "_")[[1]][1]})) %>% as.character()


# 

hom.trans.result <- 
  
  parallel::mclapply(result.list.files.genes, function(x) {
    
  treeS <- mut.files.list[grep(x, mut.files.list)]
  treeT <- wt.files.list[grep(x, wt.files.list)]
  outfile.mut_wt <- paste0(hom.trans.outfile.parent.folder,"/mut-vs-wt/", x, "_.alml")
  outfile.wt_wt <- paste0(hom.trans.outfile.parent.folder,"/wt-vs-wt/", x, "_.alml")

  
  if(Cal==T){
    
    ggvita::DELTA(
      treeS = treeS,
      treeT = treeT,
      cost = cost,
      outfile = outfile.mut_wt,
      method = "l",
      max_target = 1000,
      test = 10,
      prune = prune,
      DELTA.address = "~/app/DELTA"
    )
    
    ggvita::DELTA(
      treeS = treeT,
      treeT = treeT,
      cost = cost,
      outfile = outfile.wt_wt,
      method = "l",
      max_target = 1000,
      test = 10,
      prune= prune,
      DELTA.address = "~/app/DELTA"
    )  
    
  }
 

  mut_wt <-
    ggvita::readal(
      fileS = treeS,
      fileT = treeT,
      outfile = outfile.mut_wt,
      cost=cost,
      method = "l"
      )

  wt_wt <-
    ggvita::readal(
      fileS = treeT,
      fileT = treeT,
      outfile = outfile.wt_wt,
      cost=cost,
      method = "l"
    )

  return(list(mut_wt=mut_wt,
       wt_wt=wt_wt))

},
mc.cores = 20,
mc.silent = T)




names(hom.trans.result) <- as.character(result.list.files.genes)



hom.trans.score.list.wt_wt <- lapply(names(hom.trans.result),function(x){
  
  r <- hom.trans.result[[x]][["wt_wt"]]
  dt <- data.frame(Gene=rep(x,length(r)))
  dt$RootS <- sapply(r,function(i)i$RootS)
  dt$RootT <- sapply(r,function(i)i$RootT)
  dt$Score<- sapply(r,function(i)i$Score)
  dt$Score_order<- sapply(r,function(i)i$score_order)
  
  dt
  
  
})

names(hom.trans.score.list.wt_wt) <- result.list.files.genes


hom.trans.score.list.mut_wt <- lapply(names(hom.trans.result),function(x){
  
  r <- hom.trans.result[[x]][["mut_wt"]]
  dt <- data.frame(Gene=rep(x,length(r)))
  dt$RootS <- sapply(r,function(i)i$RootS)
  dt$RootT <- sapply(r,function(i)i$RootT)
  dt$Score<- sapply(r,function(i)i$Score)
  dt$Score_order<- sapply(r,function(i)i$score_order)
  dt$Prune.Number <- sapply(r,function(i){
    
    S.pr <- sapply(i$PruneS,function(x){
      attr(r,"params")$fileS$Lineage %>% 
        startsWith(.,x) %>% 
        sum()
      
    })
    
    T.pr <- sapply(i$PruneT,function(x){
      attr(r,"params")$fileT$Lineage %>% 
        startsWith(.,x) %>% 
        sum()
      
    })

    S.pr <- ifelse(is.na(S.pr),0,S.pr)
    T.pr <- ifelse(is.na(T.pr),0,T.pr)
    
    sum(S.pr*2-1,T.pr*2-1)
    
    })
  dt
  
  
})


names(hom.trans.score.list.mut_wt) <- result.list.files.genes




hom.trans.score.list.mut_wt.exclude <- list()


for(x in as.character(result.list.files.genes)){
  
  cat(x,"\n")
  wt_wt.df <- 
    hom.trans.score.list.wt_wt[[x]] %>% 
    #dplyr::filter(RootS!=RootT) %>% 
    dplyr::arrange(as.numeric(Score_order)) %>% 
    dplyr::mutate(Score_order=1:nrow(.))#%>% filter(Score_order<=0)
  
  mut_wt.df <- 
    hom.trans.score.list.mut_wt[[x]] %>% 
    #dplyr::filter(RootS!=RootT) %>% 
    dplyr::arrange(as.numeric(Score_order)) %>% 
    dplyr::mutate(Score_order=1:nrow(.)) %>% 
    dplyr::mutate(RootS.p = RootS %>% substr(0,nchar(.)-1),
                  RootT.p = RootT %>% substr(0,nchar(.)-1)) #%>% filter(RootS.p!=RootT.p)
  
  names(wt_wt.df)[c(4,5)] <- c("Score.wt_wt","Score_order.wt_wt")
  
  names(mut_wt.df)[c(4,5)] <- c("Score.mut_wt","Score_order.mut_wt")
  
  m.df <- merge(mut_wt.df,wt_wt.df,all = T,by=c("Gene","RootS","RootT"))
  
  m.df <- 
    m.df %>% 
    dplyr::mutate(Score_order.wt_wt=ifelse(is.na(Score_order.wt_wt),1001,Score_order.wt_wt))
  
  m.df <- 
    m.df %>% 
    dplyr::filter(is.na(Score_order.mut_wt)==F) %>% 
    dplyr::mutate(Score.order.elevation = Score_order.wt_wt-Score_order.mut_wt) %>% 
    dplyr::arrange(-as.numeric(Score.order.elevation)) %>% 
    dplyr::mutate(Elevation.order=1:nrow(.)) 
  # %>% filter(Elevation.order<=250)# %>% filter(New.Score_order<=100)
  
  
  hom.trans.score.list.mut_wt.exclude[[x]] <- m.df
  
  
}

hom.trans.score.list.mut_wt.exclude <- 
  parallel::mclapply(as.character(result.list.files.genes),function(x){
  
  wt_wt.df <- 
    hom.trans.score.list.wt_wt[[x]] %>% 
    dplyr::filter(RootS!=RootT) %>% 
    dplyr::arrange(as.numeric(Score_order)) %>% 
    dplyr::mutate(Score_order=1:nrow(.))#%>% filter(Score_order<=0)
  
  mut_wt.df <- 
    hom.trans.score.list.mut_wt[[x]] %>% 
    dplyr::filter(RootS!=RootT) %>% 
    dplyr::arrange(as.numeric(Score_order)) %>% 
    dplyr::mutate(Score_order=1:nrow(.)) %>% 
    dplyr::mutate(RootS.p = RootS %>% substr(0,nchar(.)-1),
           RootT.p = RootT %>% substr(0,nchar(.)-1)) #%>% filter(RootS.p!=RootT.p)
  
  names(wt_wt.df)[c(4,5)] <- c("Score.wt_wt","Score_order.wt_wt")
  
  names(mut_wt.df)[c(4,5)] <- c("Score.mut_wt","Score_order.mut_wt")
  
  m.df <- merge(mut_wt.df,wt_wt.df,all = T,by=c("Gene","RootS","RootT"))
  
  m.df <- 
    m.df %>% 
    dplyr::mutate(Score_order.wt_wt=ifelse(is.na(Score_order.wt_wt),1001,Score_order.wt_wt))
  
  m.df <- 
    m.df %>% 
    dplyr::filter(is.na(Score_order.mut_wt)==F) %>% 
    dplyr::mutate(Score.order.elevation = Score_order.wt_wt-Score_order.mut_wt) %>% 
    dplyr::arrange(-as.numeric(Score.order.elevation)) %>% 
    dplyr::mutate(Elevation.order=1:nrow(.)) 
  # %>% filter(Elevation.order<=250)# %>% filter(New.Score_order<=100)
      
  
  m.df
  
},mc.cores=30)


hom.trans.score.list.mut_wt.exclude <- Reduce(rbind,hom.trans.score.list.mut_wt.exclude)


hom.trans.Score_order <- hom.trans.score.list.mut_wt.exclude


hom.trans.digital <- read.table("./data/cell_to_cell_fateChange_data/new_cell_to_cell.txt",header = T,stringsAsFactors = F,colClasses = "character")
hom.trans.digital$Gene <- hom.trans.digital$Gene %>% toupper()
hom.trans.digital.uni <- hom.trans.digital[,c("from_cell_LN","to_cell_LN")] %>% unique()
hom.trans.digital.uni$Class <- c(1:nrow(hom.trans.digital.uni))
hom.trans.sister <- read.csv("~/2017-2018/tree/transform ppp to 111/symmetry-sisters.csv",header = T,colClasses = "character")


hom.trans.Score_order$Class.digital <- parallel::mclapply(1:nrow(hom.trans.Score_order),function(x){
  
  r <- hom.trans.Score_order[x,]
  
  l <- filter(hom.trans.digital,Gene==r$Gene,from_cell_LN==as.character(r$RootS),to_cell_LN==as.character(r$RootT))
  
  ifelse(nrow(l)==1,"Yes","No")

  },mc.cores = 200) %>% unlist()


hom.trans.Score_order$Class.sister <- 
  
  parallel::mclapply(1:nrow(hom.trans.Score_order),function(x){
  
  r <- hom.trans.Score_order[x,]
  
  l.1 <- filter(hom.trans.sister,Sister1_binary==as.character(r$RootS),Sister2_binary==as.character(r$RootT)) %>% nrow()
  
  l.2 <- filter(hom.trans.sister,Sister2_binary==as.character(r$RootS),Sister1_binary==as.character(r$RootT)) %>% nrow()
  
  ifelse((l.1+l.2)>=1,"Yes","No")
  
},mc.cores = 200) %>% unlist()




hom.trans.Score_order <- 
hom.trans.Score_order %>%
  mutate(Fate.Path=paste0(RootS,"-to-",RootT))


hom.trans.Score_order.FC <- 
hom.trans.Score_order %>%
  filter(#Fate.Path %in% unique(hom.trans.digital$Fate.Path),
         RootS!="Root",
         RootT!="Root",
         RootS!="RootT",
         Prune.Number<5,
         Score.order.elevation>10,
         Class.digital=="No") %>% 
  arrange(Prune.Number,-Score.order.elevation) 



hom.trans.Score_order.FC$IsFound <- 
  sapply(1:nrow(hom.trans.Score_order.FC),function(x){
    r <- hom.trans.Score_order.FC[x,]
    Genes <- Fate_change.s3[Fate_change.s3$Alternative.Path.BIN==r$Fate.Path,4]
    r$Gene %in% toupper(strsplit(as.character(Genes),split=",")[[1]])
    
    
  })



pdf("./FC-common-pr40-selected103.pdf",onefile = T,paper = "a4")
pp <- 
  parallel::mclapply(1:nrow(hom.trans.Score_order.FC),function(i){
  
  My.Fig(hom.trans.Score_order.FC[i,]$Gene,hom.trans.Score_order.FC[i,]$Score_order.mut_wt)
  
},mc.cores=20)
pp
dev.off()


# show a example 

My.Fig(gene = "C01A2.5",score.order = 1)

setDT(hom.trans.Score_order)

RootS_LN <- mclapply(hom.trans.Score_order$RootS,Bin_to_LN,mc.cores = 20)
RootT_LN <- mclapply(hom.trans.Score_order$RootT,Bin_to_LN,mc.cores = 20)

hom.trans.Score_order$RootS_LN <- RootS_LN %>% list.mapv(.[1])
hom.trans.Score_order$RootT_LN <- RootT_LN %>% list.mapv(.[1])

hom.trans.Score_order$Fate.Path.LN <- paste0(hom.trans.Score_order$RootS_LN,"-to-",hom.trans.Score_order$RootT_LN)



