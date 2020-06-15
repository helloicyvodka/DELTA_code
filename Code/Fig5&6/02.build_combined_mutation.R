library(plyr)
library(dplyr)
library(data.table)
library(parallel)
library(rlist)


mut <-list()

mut$file.list <- list.files("./data/data_from_digital/Cell_Lineage_pertution/")

mut$file.list <- mut$file.list[wt$file.list!="readme.txt"]

mut$HYP <- mut$file.list [grepl("HYP",mut$file.list )]
mut$NEU <- mut$file.list [grepl("NEU",mut$file.list )]
mut$PHA <- mut$file.list [grepl("PHA",mut$file.list )]

mut.genes <- unlist( sapply(mut$file.list,function(x){
  
  y <- strsplit(x,split="_")[[1]][1]
  
  as.character(y)
  
}) )


mut.genes <- unique(mut.genes)

mut$gene.file.list <- lapply(mut.genes,function(x)mut$file.list[grepl(x,mut$file.list)])

names(mut$gene.file.list) <- mut.genes



mut.df <- data.frame(gene=mut.genes,stringsAsFactors = F)


for(i in Tis.3){
  
  mut.df[[i]] <- sapply(mut.genes,function(x){
    
    
    any(grepl(i,mut$gene.file.list[[x]])==T)
    
    
  })
  
}




mut.df$num.marker <- sapply(1:nrow(mut.df),function(x){
  
  r <-  mut.df[x,][2:4]
  
  length(r[r==T])
    
    
  })





gene.having.3.marker <- mut.df[mut.df$num.marker==3,"gene"]





#g <- "C08B11.3"

get.mut.expr.df <- function(g){
  

  
  gene <-list()
  gene$file.list <- mut$gene.file.list[[g]]
  gene$HYP <- gene$file.list [grepl("HYP",gene$file.list )]
  gene$NEU <- gene$file.list [grepl("NEU",gene$file.list )]
  gene$PHA <- gene$file.list [grepl("PHA",gene$file.list )]
  
  
  gene.dt <- lapply(c("HYP", "NEU", "PHA"), function(m) {
    d <- lapply(gene[[m]], function(x) as.data.frame(fread(paste0("./data/data_from_digital/Cell_Lineage_pertution/",x),stringsAsFactors = F)))
    names(d) <- gene[[m]]
    for(i in 1:length(d)){
      names(d[[i]])[names(d[[i]])=="binary_call"] <- paste0(m,"_",i)
      
    }
    d
  })
  
  names(gene.dt) <- c("HYP","NEU","PHA")
  
  Get.marker.expr <- function(m){
    
    The.full_join <- function(x,y){full_join(as.data.frame(x),as.data.frame(y)[,c(1,3)],by="cell_name")}
    
    dt.1 <- Reduce(The.full_join,gene.dt[[m]])[,c(-2,-4)]
    
    dt.1[[m]] <- sapply(1:nrow(dt.1),function(x){
      
      r <- dt.1[x,-1] %>% as.character()
      r <- r[r!="NA"] %>% na.omit()
      ifelse((r[r=="1"] %>% length())>= 
             (r[r=="0"] %>% length()),
             1,
             0)
      
      
      
    })
    
    dt.1[,c(1,ncol(dt.1))]
  }
  
  
  
  
Get.com.expr <- function(l){
    
    the.List.1 <- lapply(l, Get.marker.expr)
    Reduce(function(x,y)full_join(x,y,by="cell_name"),the.List.1)
    
  }
  
  

  gene.com.expr <- Get.com.expr(list("HYP","NEU","PHA" ))
  
  gene.com.expr <- gene.com.expr %>% na.omit() 
  
  gene.com.expr$Lineage <- gene.com.expr$cell_name %>% LN_to_Bin()
  
  # HYP=1, NEU=2, PHA=4
  
  
  gene.com.expr$Class.num<- sapply(1:nrow(gene.com.expr),function(x){
    
    n <- c(1,2,4)
    r <- gene.com.expr[x,2:4]
    
    y <- n[r==1] %>% as.integer()
    sum(y)
    
  })
 
  
  The.Lineage <- gene.com.expr$Lineage 
  The.Lineage.2 <- c("Root",
                     "0","1",
                     "00","01","10","11",
                     "000","001","010","011","100","101","110","111")
  
  
  The.addNode <- setdiff(The.Lineage.2,The.Lineage)
  The.addNode.df <- data.frame(Lineage=The.addNode,stringsAsFactors = F)
  gene.com.expr <- full_join(The.addNode.df,gene.com.expr,by="Lineage")

  
  
  
  repeat{
  if(
    
    gene.com.expr$Lineage %>% as.character() %>%  Find.addNode() %>% length()!=0
    
  ){
    
    gene.com.expr <- gene.com.expr %>% Na.omit.expr.dt()
    
  }else{
    break
    }
  
  }
  

  return(gene.com.expr)
  
}




mut.df.list <- list() 
mut.full.df.list <- list() 
for(g in gene.having.3.marker ){
  cat(g,"\n")
  mut.df.list[[g]]<-get.mut.expr.df(g)
  mut.full.df.list[[g]]<-Na.omit.expr.dt(mut.df.list[[g]])
}


# Gene N
# 1: ATX-2 2
# 2: CUL-3 2
# 3: MOG-4 2
# 4: PAR-5 1
# 5: PKC-3 1
# 6: RAN-4 2





