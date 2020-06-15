library(plyr)
library(dplyr)
library(data.table)
library(parallel)
library(rlist)



wt <-list()
wt$file.list <- list.files("./data/data_from_digital/WT_20180604/",full.names = T)
wt$file.list <- wt$file.list[wt$file.list!="readme.txt"]
wt$HYP <- wt$file.list [grepl("RW10348",wt$file.list )]
wt$NEU <- wt$file.list [grepl("RW10434",wt$file.list )]
wt$PHA <- wt$file.list [grepl("RW10425",wt$file.list )]


wt.dt <- lapply(c("HYP", "NEU", "PHA"), function(m) {
  d <- lapply(wt[[m]], function(x)
    fread(x))
  names(d) <- wt[[m]]
  d
})
names(wt.dt) <- c("HYP","NEU","PHA")

# # com = common
# 
# 
# com.cell.HYP <- Reduce(union,lapply(wt.dt$HYP,function(x)x$cell_name))
# com.cell.NEU <- Reduce(union,lapply(wt.dt$NEU,function(x)x$cell_name))
# com.cell.PHA <- Reduce(union,lapply(wt.dt$PHA,function(x)x$cell_name))
# 
# com.cell <- Reduce(union,list(com.cell.HYP,com.cell.NEU ,com.cell.PHA ))



for(t in c("HYP","NEU","PHA")){
  for(i in 1:10){
    names(wt.dt[[t]][[i]])[names(wt.dt[[t]][[i]])=="binary_call"] <- paste0(t,"_",i)
    
  }
}


com.dt <- list()

for(t in c("HYP","NEU","PHA")){
  
  dt <- Reduce(function(x,y){full_join(as.data.frame(x),as.data.frame(y)[,c(1,3)],by=c("cell_name"))},wt.dt[[t]])
  
  dt <- dt[,!names(dt) %in% c("average_expression","time_point_expression") ]
  
  #dt <- dt %>% mutate(Lineage=LN_to_Bin(cell_name))
  
  #dt <- data.frame(Lineage=LN_to_Bin(dt$cell_name),dt[,1:ncol(dt)])
  
  # dt <- dt[ dt$cell_name %in% com.cell,]
  
  # assign(paste0("com.dt.",t),dt)
  
  com.dt[[t]] <- dt
  
  write.table(dt,file=paste0("./data/data_from_digital/WT_merge/20180604_",t,".txt"),quote = F,row.names = F)
  
  rm(dt)
}


for(t in c("HYP","NEU","PHA")){
  
  dt <- com.dt[[t]]
  
  dt[[paste0("all.",t)]] <- sapply(1:nrow(dt),function(x){
    
    r <- dt[x,2:ncol(dt)]
    
    if(length(r[r==1])>=length(r[r==0])){
      
      y<-1
    }else{
      
      y<-0
    }
    
    y
    
  })
  
  com.dt[[t]] <-  dt
  
  rm(dt)
  
}

com.dt[[1]] <- com.dt[[1]][,c(1,ncol(com.dt[[1]]))]

com.expr <-  Reduce(function(x,y){full_join(as.data.frame(x),as.data.frame(y)[c(1,ncol(y))],by=c("cell_name"))},com.dt)



Tis.3 <- c("HYP","NEU","PHA")

# HYP=1, NEU=2, PHA=4
com.expr$Class.num<- sapply(1:nrow(com.expr),function(x){
  
  n <- c(1,2,4)
  r <- com.expr[x,2:4]
  
  

  y <- n[r==1]
  sum(y)
  
})


com.expr$Lineage <- com.expr$cell_name %>% LN_to_Bin()

Na.omit.expr.dt <- function(com.expr){
  
  com.expr.tips <- com.expr %>% filter(Lineage %in% (Lineage %>% Find.tips())) %>% na.omit() %>% `$`("Lineage")
  
  com.expr.tips.na <- setdiff(com.expr %>% filter(Lineage %in% (Lineage %>% Find.tips()))  %>% select(Lineage) %>% `$`("Lineage"),
                              com.expr.tips)
  
  com.expr.2 <- com.expr  %>% filter(!Lineage %in% com.expr.tips.na)
  
  if(com.expr.2$Lineage %>% Find.addNode() %>% length()!=0)stop("Error: Tree is not complete!")
  
  return(com.expr.2)
  
}

com.expr.2 <- Na.omit.expr.dt(com.expr)




#write.table(com.expr.2,file=paste0("./data/data_from_digital/WT_merge/","20180612_HYP_NEU_PHA",".txt"),quote = F,row.names = F)
