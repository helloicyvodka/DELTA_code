wt.2$Lineage.parent <- sapply(wt.2$Lineage,get_parent)

wt.2 <- wt.2[!wt.2$Lineage %in% wt.2$Lineage.parent,]

wt.3 <- wt.2[,c(1,2,3)]

write.table(wt.3,file="../WT_merge/wt_for_DELTA.alm",quote = F,row.names = F)

rm_tail_blank_line("../WT_merge/wt_for_DELTA.alm")



# write full mutated tree 


mclapply(names(mut.df.list),function(x){
  
  df <- mut.df.list[[x]]
  
  df$Lineage <- LN_to_Bin(df$cell_name)
  
  df.2 <- data.frame(Lineage=df$Lineage,Name=df$cell_name,Class=df$Class.num)
  
  file.name <-  paste0("../mutated_for_DELTA/mutated_full_for_DELTA/",x,"_for_DELTA.txt") 
  
  write.table(wt.2,file=file.name,quote = F,row.names = F)
  
  rm_tail_blank_line(file.name)
  
},mc.cores = 20)






# write mutated tree with only tips


mclapply(names(mut.df.list),function(x){
  
  df <- mut.df.list[[x]]
  
  df$Lineage <- LN_to_Bin(df$cell_name)
  
  df$Lineage.parent <- sapply(df$Lineage,get_parent)
  
  df <- df[ ! df$Lineage %in% df$Lineage.parent,]
  
  df.2 <- data.frame(Lineage=df$Lineage,Name=df$cell_name,Class=df$Class.num)
  
  file.name <-  paste0("../mutated_for_DELTA/mutated_for_DELTA/",x,"_for_DELTA.alm") 
  
  write.table(df.2,file=file.name,quote = F,row.names = F)
  
  rm_tail_blank_line(file.name)
  
},mc.cores = 20)




# write wt cut tree




















##########################################################

### Date : 2018-07-02

for(x in names(mut.full.df.list)){
  
  cat(x,"\n")
  
  mut.1 <- mut.full.df.list[[x]]
  wt <- com.expr.2  
  
  Com.Lineage <- intersect(wt$Lineage,mut$Lineage)
  
  mut <- mut %>% filter(Lineage %in% Com.Lineage) %>% filter(Lineage %in% (Lineage %>% Find.tips()))
  wt  <- wt  %>% filter(Lineage %in% Com.Lineage) %>% filter(Lineage %in% (Lineage %>% Find.tips()))
  
  mut.2 <-data.frame(Lineage=mut$Lineage,Name=mut$cell_name,Class=mut$Class.num,stringsAsFactors = F)
  wt.2 <- data.frame(Lineage=wt$Lineage,Name=wt$cell_name,Class=wt$Class.num,stringsAsFactors = F)
  
  
  wt.file.name <-  paste0("./data2/wt-tree-tips/",x,"_.alm")
  mut.file.name <-  paste0("./data2/mut-tree-tips/",x,"_.alm")
  
  write.alm(wt.2,file=wt.file.name)
  write.alm(mut.2,file=mut.file.name)
  
  
}


BioGRID.3$DELTA.score <-
  parallel::mclapply(1:nrow(BioGRID.3),function(i){
  
  r <- BioGRID.3[i,]
  
  output <- paste0("./data3/gene-network/",i)
  
  dir.create(output)
  
  
  
  mut.A <- mut.full.df.list[[r$Index.A.2]]
  mut.B <- mut.full.df.list[[r$Index.B.2]]
  
  if(is.null(mut.B)|is.null(mut.A)){
    
    return(NA)
    
  }else{
    
    AB.inter <- intersect(mut.A$Lineage,mut.B$Lineage) 
    
    AB.inter <- c(c("Root","0","1","00","01","11","10"),AB.inter) %>% unique()
    
    if(is.null(ggvita::Find.addNode(AB.inter))){
      
      AB.tips <- ggvita::Find.tips(AB.inter)
      
      mut.A.2 <-data.frame(Lineage=mut.A$Lineage,Name=mut.A$cell_name,Class=mut.A$Class.num,stringsAsFactors = F)
      mut.B.2 <-data.frame(Lineage=mut.B$Lineage,Name=mut.B$cell_name,Class=mut.B$Class.num,stringsAsFactors = F)
      
      ggvita::write.alm(mut.A.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$Index.A.2,".alm"))
      ggvita::write.alm(mut.B.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$Index.B.2,".alm"))
      
      DELTA(
        treeS = paste0(output,"/",r$Index.A.2,".alm"),
        treeT = paste0(output,"/",r$Index.B.2,".alm"),
        cost = cost,
        outfile = paste0(output,"/","My.almg"),
        method = "g",
        #max_target = 1000,
        test = 10
      )
      
      t <- ggvita::readal.almg(paste0(output,"/","My.almg"))
      
      return(t$Score)
    }
    
  }
  
  
  
},mc.cores=20) %>% unlist()

BioGRID.4 <- BioGRID.3 %>% filter(!is.na(DELTA.score))


sapply(list.files("./data3/gene-network/",full.names = T),function(x)list.files(x,full.names = T) %>% length)



combn(names(mut.full.df.list),2,simplify = T) %>% t()

BioGRID.6 <- data.frame(combn(names(mut.full.df.list),2,simplify = T) %>% t(),stringsAsFactors = F)
colnames(BioGRID.6) <- c("Index.A.2","Index.B.2")





BioGRID  <- 
  read.delim("./data/Biogrid/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.4.163.tab2.txt",sep="\t",header = T,stringsAsFactors = F)

# BioGRID.7 <- BioGRID %>% filter(Experimental.System!="Two-hybrid"&Experimental.System.Type=="physical")

BioGRID.7<- BioGRID %>% mutate(Official.Symbol.Interactor.A = Official.Symbol.Interactor.A %>% toupper(),
                               Official.Symbol.Interactor.B = Official.Symbol.Interactor.B %>% toupper())

BioGRID.7$Index.A <- sapply(1:nrow(BioGRID.7),function(x){
  
  r <- BioGRID.7[x,]
  My.all <- list(r$Systematic.Name.Interactor.A,r$Official.Symbol.Interactor.A,r$Synonyms.Interactor.A) %>% unlist() %>% as.character()
  My.all <- sapply(My.all,function(i)strsplit(i,split='\\|')) %>% unlist(recursive = T) %>% as.character()%>% gsub("CELE_","",x = .) 
  My.all[My.all!="-"] %>% unique() %>% toupper()
})

BioGRID.7$Index.B <- sapply(1:nrow(BioGRID.7),function(x){
  
  r <- BioGRID.2[x,]
  My.all <- list(r$Systematic.Name.Interactor.B,r$Official.Symbol.Interactor.B,r$Synonyms.Interactor.B) %>% unlist() %>% as.character()
  My.all <- sapply(My.all,function(i)strsplit(i,split='\\|')) %>% unlist(recursive = T) %>% as.character()%>% gsub("CELE_","",x = .) 
  My.all[My.all!="-"] %>% unique()%>% toupper()
})




BioGRID.7$Index.A.2 <- sapply(BioGRID.7$Index.A,function(x){
  
  x[x %in% mut.genes]
  
})

BioGRID.7$Index.B.2 <- sapply(BioGRID.7$Index.B,function(x){
  
  x[x %in% mut.genes]
  
})



My.test.1 <- setdiff(BioGRID.6[,c(1,2)] %>% rlist::list.parse(),BioGRID.7[,c("Index.A.2","Index.B.2")] %>% rlist::list.parse())
My.test.2 <- setdiff(BioGRID.6[,c(1,2)] %>% rlist::list.parse(),BioGRID.7[,c("Index.B.2","Index.A.2")] %>% rlist::list.parse())

BioGRID.8 <- union(My.test.1,My.test.2) %>% rlist::list.stack()

BioGRID.9 <- BioGRID.8[sample(1:nrow(BioGRID.8),30),]

rownames(BioGRID.9) <- 1:30

BioGRID.9$DELTA.score <-
  parallel::mclapply(1:nrow(BioGRID.9),function(i){
    
    r <- BioGRID.9[i,]
    
    output <- paste0("./data3/gene-network2/",i)
    
    dir.create(output)
    
    
    
    mut.A <- mut.full.df.list[[r$Index.A.2]]
    mut.B <- mut.full.df.list[[r$Index.B.2]]
    
    if(is.null(mut.B)|is.null(mut.A)){
      
      return(NA)
      
    }else{
      
      AB.inter <- intersect(mut.A$Lineage,mut.B$Lineage) 
      
      AB.inter <- c(c("Root","0","1","00","01","11","10"),AB.inter) %>% unique()
      
      if(is.null(ggvita::Find.addNode(AB.inter))){
        
        AB.tips <- ggvita::Find.tips(AB.inter)
        
        mut.A.2 <-data.frame(Lineage=mut.A$Lineage,Name=mut.A$cell_name,Class=mut.A$Class.num,stringsAsFactors = F)
        mut.B.2 <-data.frame(Lineage=mut.B$Lineage,Name=mut.B$cell_name,Class=mut.B$Class.num,stringsAsFactors = F)
        
        ggvita::write.alm(mut.A.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$Index.A.2,".alm"))
        ggvita::write.alm(mut.B.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$Index.B.2,".alm"))
        
        DELTA(
          treeS = paste0(output,"/",r$Index.A.2,".alm"),
          treeT = paste0(output,"/",r$Index.B.2,".alm"),
          cost = cost,
          outfile = paste0(output,"/","My.almg"),
          method = "g",
          #max_target = 1000,
          test = 10
        )
        
        t <- ggvita::readal.almg(paste0(output,"/","My.almg"))
        
        return(t$Score)
      }
      
    }
    
    
    
  },mc.cores=20) %>% unlist()


library(ggplot2)
library(cowplot)


data.BioGRID <- rbind(BioGRID.4[,c("Index.A.2","Index.B.2","DELTA.score")] %>% mutate(group="BioGRID"),
                      BioGRID.9[,c("Index.A.2","Index.B.2","DELTA.score")] %>% mutate(group="non-BioGRID"))

data.BioGRID.2  <- data.BioGRID %>% group_by(group) %>% summarise(Score.mean=mean(as.numeric(DELTA.score)),
                                                                  Score.se =std(as.numeric(DELTA.score)))

p.B <-
  ggplot(data.BioGRID.2)+
  geom_errorbar(aes(x=group,ymax=Score.mean+Score.se,ymin=Score.mean-Score.se),width=0.1)+
  geom_point(aes(x=group,y=Score.mean),shape=18,size=3)+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("DELTA score")+
  xlab(NULL)
  
  
Test <- wilcox.test(BioGRID.4$DELTA.score %>% as.numeric(),BioGRID.9$DELTA.score %>% as.numeric(),paired = F)

label <- substitute(paste("Wilcox test  ", " P = ", pvalue),
                    list(pvalue = signif(Test$p.value, 2)))


ggdraw(p.B) + draw_label(label, .6, .9,size=10) 













