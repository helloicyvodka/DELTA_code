library(ggplot2)
library(magrittr)
library(dplyr)
setwd("~/2017-2018/tree/tree_fate_change")


second_name <- xlsx::read.xlsx("~/2017-2018/tree/tree_fate_change/revised_2020/1-s2.0-S1534580715004876-mmc4.xlsx",sheetIndex = 3)
second_name <- second_name[!is.na(second_name$Other_Names),]
second_name_dict <- second_name$Gene_Name %>% toupper()
names(second_name_dict) <- second_name$Other_Names %>% toupper()
second_name_dict <- c(second_name_dict,"SKR-1/2"="SKR-2","OMA-1/2"="OMA-1")
To_second_name <- function(x){
  x <- as.character(x)
  if(x %in% names(second_name_dict)) x <- second_name_dict[x] %>% as.character()
  x
}


hom.trans.digital <- read.table("~/2017-2018/tree/tree_fate_change/data/cell_to_cell_fateChange_data/new_cell_to_cell.txt",header = T,stringsAsFactors = F,colClasses = "character")
hom.trans.digital$Gene <- hom.trans.digital$Gene %>% toupper()
hom.trans.digital$Gene <- hom.trans.digital$Gene %>% sapply(To_second_name)

#hom.trans.digital$Gene

mut.files.list <- list.files("../tree_fate_change/data2/mut-tree-tips/",full.names = T)
wt.files.list <- list.files("../tree_fate_change/data2/wt-tree-tips/",full.names = T)
cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"


out.folder.all <- "~/2017-2018/tree/tree_fate_change/data2/hom.trans.pvalue.20200420/"
dir.create(out.folder.all)

hom.trans.digital$pvalue <- 
  parallel::mclapply(1:nrow(hom.trans.digital),function(x){
    
    r <- hom.trans.digital[x,]
    
    if(length(mut.files.list[grep(r$Gene, mut.files.list)])==0|length(wt.files.list[grep(r$Gene,wt.files.list)])==0){
      return("No three tissue markers")
    }else{
        
      out.folder <- paste0(out.folder.all,as.character(x))
      outfile<- paste0(out.folder, "/global.almg")
      # 
      # 
      # dir.create(out.folder)
      # 
      # 
      # 
      # treeS <-  mut.files.list[grep(r$Gene, mut.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
      # treeS <- treeS %>% filter(Lineage %>% startsWith(r$from_cell_LN))
      # treeS$Lineage <- treeS$Lineage %>% sub(r$from_cell_LN,"",.)
      # 
      # 
      # 
      # treeT <-  wt.files.list[grep(r$Gene, wt.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
      # treeT <- treeT %>% filter(Lineage %>% startsWith(r$to_cell_LN))
      # treeT$Lineage <- treeT$Lineage %>% sub(r$to_cell_LN,"",.)
      # 
      # #if(nrow(treeS)<=4 | nrow(treeT)<=4)return(NA)
      # 
      # 
      # 
      # write.alm(treeS,paste0(out.folder,"/treeS.alm"))
      # write.alm(treeT,paste0(out.folder,"/treeT.alm"))
      # 
      # 
      # ggvita::DELTA(
      #   treeS = paste0(out.folder,"/treeS.alm"),
      #   treeT = paste0(out.folder,"/treeT.alm"),
      #   cost = cost,
      #   outfile = outfile,
      #   method = "g",
      #   test = 1000,
      #   prune = 1,
      #   all="F"
      # )


      rr <- ggvita::readal.almg(outfile)
      
      
      
      return(rr$PValue$pvalue)
      
      
      }
    
  },
  mc.cores = 20)



hom.trans.digital$Score <- 
  parallel::mclapply(1:nrow(hom.trans.digital),function(x){
    
    r <- hom.trans.digital[x,]
    
    if(length(mut.files.list[grep(r$Gene, mut.files.list)])==0|length(wt.files.list[grep(r$Gene,wt.files.list)])==0){
      return("No three tissue markers")
    }else{
      
      out.folder <- paste0(out.folder.all,as.character(x))
      outfile<- paste0(out.folder, "/global.almg")
      # 
      # 
      # dir.create(out.folder)
      # 
      # 
      # 
      # treeS <-  mut.files.list[grep(r$Gene, mut.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
      # treeS <- treeS %>% filter(Lineage %>% startsWith(r$from_cell_LN))
      # treeS$Lineage <- treeS$Lineage %>% sub(r$from_cell_LN,"",.)
      # 
      # 
      # 
      # treeT <-  wt.files.list[grep(r$Gene, wt.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
      # treeT <- treeT %>% filter(Lineage %>% startsWith(r$to_cell_LN))
      # treeT$Lineage <- treeT$Lineage %>% sub(r$to_cell_LN,"",.)
      # 
      # #if(nrow(treeS)<=4 | nrow(treeT)<=4)return(NA)
      # 
      # 
      # 
      # write.alm(treeS,paste0(out.folder,"/treeS.alm"))
      # write.alm(treeT,paste0(out.folder,"/treeT.alm"))
      # 
      # 
      # ggvita::DELTA(
      #   treeS = paste0(out.folder,"/treeS.alm"),
      #   treeT = paste0(out.folder,"/treeT.alm"),
      #   cost = cost,
      #   outfile = outfile,
      #   method = "g",
      #   test = 1000,
      #   prune = 1,
      #   all="F"
      # )
      
      
      rr <- ggvita::readal.almg(outfile)
      
      
      
      return(rr$Score)
      
      
    }
    
  },
  mc.cores = 20)




#data.table::fwrite(hom.trans.digital,"./revised_2020/hom_trans_digital.csv")


rr.2 <- 
  parallel::mclapply(1:nrow(hom.trans.digital), function(x) {
    
    r <- hom.trans.digital[x,]
    out.folder <- paste0("./data2/hom.trans.pvalue.20200420/",as.character(x))
    outfile<- paste0(out.folder, "/global.almg")
    
    
    if(length(mut.files.list[grep(r$Gene, mut.files.list)])==0|length(wt.files.list[grep(r$Gene,wt.files.list)])==0){return(NA)}
    
    rr <- ggvita::readal.almg(outfile)
    
    
    
    return(rr)
    
  },
  mc.cores = 10)


hom.trans.digital$TipSize <-
  rr.2 %>% sapply(function(x){
    if(is.na(x))return(NA)
    else{return(c(x$MatchS,x$RootS) %>% ggvita::Find.tips() %>% length())}
    
  })


# ggplot(hom.trans.digital %>% na.omit())+
#   geom_histogram(aes(-log(as.numeric(pvalue),base=10)),binwidth = 1,fill="#8EE5EE")+
#   scale_x_continuous(breaks =seq(0,100,10),expand = c(0,0),limits = c(0,90))+
#   scale_y_continuous(breaks =seq(0,10,1),expand = c(0.005,0),limits = c(0,11))+
#   labs(title="Homeotic transformation p-value test",x="-log10(P-value)",y="Density")+
#   theme_bw()+
#   geom_text(label="N = 141",x=50,y=10)+
#   geom_segment(x = -log10(0.05),xend = -log10(0.05),y=-0.5,yend=11,size=0.5,color="red")+
#   theme(axis.text.x = element_text(size=5),
#         plot.title  = element_text(hjust = 0.5),
#         plot.background = element_rect(fill = NULL,color = NULL),panel.grid = element_blank())
# 



# 
# rr <- list()
# 
# 
# for(x in 1:nrow(hom.trans.digital)){
# 
#   cat(x,"\n")
# 
#   r <- hom.trans.digital[x,]
# 
#   out.folder <- paste0("./data2/hom.trans.pvalue.2/",as.character(x))
# 
#   if(length(mut.files.list[grep(r$Gene, mut.files.list)])==0|length(wt.files.list[grep(r$Gene,wt.files.list)])==0){rr[[x]]<-NA;next}
# 
# 
#   dir.create(out.folder)
# 
# 
# 
#   treeS <-  mut.files.list[grep(r$Gene, mut.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
#   treeT <-  wt.files.list[grep(r$Gene, wt.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
# 
# 
# 
#   treeS <- treeS %>% filter(Lineage %>% startsWith(r$from_cell_LN))
#   treeS$Lineage <- treeS$Lineage %>% sub(r$from_cell_LN,"",.)
# 
# 
#   treeT <- treeT %>% filter(Lineage %>% startsWith(r$to_cell_LN))
#   treeT$Lineage <- treeT$Lineage %>% sub(r$to_cell_LN,"",.)
# 
#   
#   
#   # c(treeS$Lineage %>% ggvita::Find.tree.from.tips() %>% ggvita::Find.addNode(),
#   #   treeT$Lineage %>% ggvita::Find.tree.from.tips() %>% ggvita::Find.addNode()) %>% 
#   #   print()
#   
#   
#   #write.alm(treeS,paste0(out.folder,"/treeS.alm"))
#   #write.alm(treeT,paste0(out.folder,"/treeT.alm"))
#   outfile<- paste0(out.folder, "/global.almg")
# 
# 
#   si <-
#     ggvita::DELTA(
#     treeS = paste0(out.folder,"/treeS.alm"),
#     treeT = paste0(out.folder,"/treeT.alm"),
#     cost = cost,
#     outfile = outfile,
#     method = "g",
#     test = 100,
#     prune = 1,
#     all="F",
#     DELTA.address = DELTA.address
#   )
#   rr[[x]] <- ggvita::readal.almg(outfile)
# 
# 
# }






# ggplot(hom.trans.digital %>% na.omit())+
# geom_histogram(aes(-log(pvalue,base=10)),binwidth = 1,fill="#8EE5EE")+
# scale_x_continuous(breaks =seq(0,100,10),expand = c(0,0),limits = c(0,90))+
# scale_y_continuous(breaks =seq(0,10,1),expand = c(0.005,0),limits = c(0,11))+
# labs(title="Homeotic transformation p-value test",x="-log10(P-value)",y="Density")+
# theme_bw()+
# geom_text(label="N = 141",x=50,y=10)+
# geom_segment(x = -log10(0.05),xend = -log10(0.05),y=-0.5,yend=11,size=0.5,color="red")+
# theme(axis.text.x = element_text(size=5),plot.title  = element_text(hjust = 0.5),plot.background = element_rect(fill = NULL,color = NULL),panel.grid = element_blank())


# hom.trans.digital.2 <- 
#   hom.trans.digital %>% 
#   na.omit() %>% 
#   mutate(pvalue.2=-log10(pvalue))  %>% 
#   arrange(pvalue.2) %>% 
#   mutate(density=(0:140)/141)
# 
# 
# ggplot(hom.trans.digital.2)+
#   geom_line(aes(x=pvalue.2,y=density),color="blue",size=1)+
#   scale_x_continuous(expand = c(0.0,0.0),limits = c(0,10.5),breaks = c(0:10))+
#   scale_y_continuous(expand = c(0.0,0.0),limits = c(0,0.4))+
#   labs(title="Homeotic transformation p-value test",x="-log10(P-value)",y="Cumulative density")+
#   geom_text(label="N = 141",x=50,y=10)+
#   geom_segment(x = -log10(0.05),xend = -log10(0.05),y=-0.5,yend=11,size=0.5,color="red")+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=10),plot.title  = element_text(hjust = 0.5),plot.background = element_rect(fill = NULL,color = NULL),panel.grid = element_blank(),panel.border = element_rect(size=1))
# 




###############   hom.trans.digital.3 hom.trans.digital.4   #################






#hom.trans.digital.3 <- hom.trans.digital.3 %>% filter(!is.na(pvalue))






# rr <- 
#   
#   parallel::mclapply(1:nrow(hom.trans.digital.3), function(x) {
#     
#     r <- hom.trans.digital.3[x,]
#     
#     if(length(mut.files.list[grep(r$Gene, mut.files.list)])==0|length(wt.files.list[grep(r$Gene,wt.files.list)])==0){return(NA)}
#     
#     treeS <-  mut.files.list[grep(r$Gene, mut.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
#     treeS <- treeS %>% filter(Lineage %>% startsWith(r$from_cell_LN))
#     treeS$Lineage <- treeS$Lineage %>% sub(r$from_cell_LN,"",.)
#     
#     treeT <-  wt.files.list[grep(r$Gene, wt.files.list)] %>% read.table(header = T,stringsAsFactors = F,colClasses = "character")
#     treeT <- treeT %>% filter(Lineage %>% startsWith(r$to_cell_LN))
#     treeT$Lineage <- treeT$Lineage %>% sub(r$to_cell_LN,"",.)
#     
#     rr <- list(Class.num.S=unique(treeS$Class) %>% length,Class.num.T=unique(treeT$Class) %>% length)
#     
#     
#     return(rr)
#     
#   },
#   mc.cores = 10)



# hom.trans.digital.3$Class.num.S <-
#   
#   rr %>% sapply(function(x){
#   
#   if(is.na(x))return(NA)
#   else{return(x$Class.num.S)}
#   
# })
# 
# 
# 
# hom.trans.digital.3$Class.num.T <-
#   
#   rr %>% sapply(function(x){
#     
#     if(is.na(x))return(NA)
#     else{return(x$Class.num.T)}
#     
#   })
# 
# hom.trans.digital.3$TipSize <-
#   rr.2 %>% sapply(function(x){
#     if(is.na(x))return(NA)
#     else{return(c(x$MatchS,x$RootS) %>% ggvita::Find.tips() %>% length())}
#     
#   })
# 
# hom.trans.digital.3$Score <-
#   rr.2 %>% sapply(function(x){
#     if(is.na(x))return(NA)
#     else{return(x$Score %>% as.numeric())}
#     
#   })


#M.rank <- nrow(hom.trans.digital)





# hom.trans.digital.4 <- 
#   hom.trans.digital.3 %>% 
#   na.omit() %>% 
#   arrange(neg.log10.pvalue) %>% 
#   mutate(density=(0:(M.rank-1))/M.rank ,rank=1:M.rank,Score.divid.TipSize=Score/TipSize)


hom.trans.digital.4 <- hom.trans.digital
hom.trans.digital.4[is.na(hom.trans.digital$pvalue),]$pvalue <- 1

hom.trans.digital.4 <- hom.trans.digital.4 %>% na.omit()

hom.trans.digital.4$pvalue <- hom.trans.digital.4$pvalue %>%  as.numeric()
hom.trans.digital.4$Score <- hom.trans.digital.4$Score %>% as.numeric()


hom.trans.digital.4 <- 
  hom.trans.digital.4 %>% 
  mutate(neg.log10.pvalue=-log10(as.numeric(pvalue)))


PPP.1 <-
ggplot(hom.trans.digital.4)+
  geom_point(aes(x=rank(TipSize,ties.method ="random"),y=Score),size=0.5,shape=21)+
  theme_bw()+
  theme(panel.background = element_blank())+
  ylab("DELTA Score")+
  xlab("Rank of tree size")+
  theme(panel.grid = element_blank())

PPP.2 <-  
ggplot(hom.trans.digital.4)+
  geom_point(aes(x=rank(TipSize,ties.method ="random" ),y=neg.log10.pvalue),size=0.5,shape=21)+
  theme_bw()+
  theme(panel.background = element_blank())+
  geom_hline(yintercept =-log10(0.05),color="red")+
  ylab("-log10(P)")+
  xlab("Rank of tree size")+
  theme(panel.grid = element_blank())


# total events = 135

cowplot::plot_grid(plotlist = list(PPP.1,PPP.2),nrow=1,labels = c("B","C"),rel_widths = c(1.8,1.7))


sum(hom.trans.digital$pvalue<=0.05,na.rm = T)/157 #87.3%
