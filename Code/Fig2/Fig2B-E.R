
rm(list = ls() )

library(plyr)
library(dplyr)
library(data.tree)
library(rlist)
library(pipeR)
library(testit)  
library(ggplot2)
library(parallel)
library(data.table)
library(cowplot)
library(ggvita)

method="l"
std <- function(x){x <- na.omit(x); sd(x)/sqrt(length(x))}

## Define a function to calculate the simulated trees
MyPlot <- function(score="common",
                   m=1,
                   u=1,
                   prune=20,
                   test.num=10,
                   Cal=F,
                   D1=0.35
){
  
  for(file in list.files("/mnt/data/home/yuanmeng/2017-2018/tree/SimulateTree/code/",full.names = T,recursive = T))source(file)
  

  #B <- c(0.05,0.10,0.15,0.20)
  B <- c(0.05,0.05,0.05,0.05)
  
  D <- c(D1,0.40,0.45,0.50)
  
  
  
  Parm.df <- data.frame(R=c(D,0,0,0,0),B=c(0,0,0,0,B),stringsAsFactors = F)
  
  Parm.df <- rbind(c(0,0),Parm.df,c(D1,0.05))
  
  Parm.df$Group <- sapply(as.character(1:nrow(Parm.df)),function(x)paste0("G",x))
  
  
  
  My.Test.df.2 <- function(test.num,R,B,parent.folder="./",Cal=Cal){
    
    outputfolder <- paste0(parent.folder,"data.R",R,"B",B,"/")    #./data.R0B0/
    
    dir.create(outputfolder)
    
    T.test <- 1:as.integer(test.num)
    
    method="l"
    
    a = 100
    N = 8
    K = 4
    Dmax = 10
    tmax = 50
    bias = as.numeric(B)
    RdropRate =as.numeric(R)
    
    
    
    if(Cal){
      
      si <- mclapply(T.test ,function(x){
        
        
        si <- F.main.2(ID =x,
                       folder = outputfolder ,
                       a = a,N=N,K=K,Dmax = Dmax,tmax = tmax,
                       bias = as.numeric(B),
                       RdropRate =as.numeric(R),
                       leafTypeMin = 6,
                       leafCountMin = 500,
                       tryMax = 100,
                       digits = 3,
                       score = score,
                       m=m,
                       u=u,
                       silent = F
        )
        
        
        
        Info  <-  paste0("ID", x,
                         "a", a,
                         "N", N,
                         "K", K,
                         "t", tmax,
                         "D", Dmax,
                         "b", bias,
                         "R", RdropRate)
        
        
        sil <- ggvita::DELTA(treeS=paste0(outputfolder,Info,"/",Info,".1.alm"),
                             treeT=paste0(outputfolder,Info,"/",Info,".2.alm"),
                             cost =paste0(outputfolder,Info,"/",Info,".tsv"),
                             outfile = paste0(outputfolder,Info,"/",Info,".alm",method),
                             method = method,
                             max_target = 10,
                             prune = prune,
                             test = 100,
                             DELTA.address = NULL)
        
        
        
      },mc.cores = 30)
      
      
    }
    
    
    Result.Df <- F.ReadResult.2(outputfolder)
    
    return(Result.Df)
    
  }
  
  
  
  
  My.Results.3 <- list()
  
  # 
  # for(i in c(1,2,9,10)){
  #   
  #   cat(i,"\n")
  #   
  #   the.line <- Parm.df[i,]
  #   
  #   Df  <-  My.Test.df.2(test.num = 10,R = the.line$R,B = the.line$B)
  #   
  #   Df$Group <- rep(the.line$Group,nrow(Df))
  #   
  #   My.Results.3[[i]] <- Df
  #   
  #   
  # }
  
  My.Results.3 <-
    parallel::mclapply(c(1,2,9,10),function(i){
      
      #cat(i,"\n")
      
      the.line <- Parm.df[i,]
      
      Df  <-  My.Test.df.2(test.num =test.num,R = the.line$R,B = the.line$B,Cal=Cal)
      
      Df$Group <- rep(the.line$Group,nrow(Df))
      
      Df
    },mc.cores = 4)
  
  
  
  
  
  My.Results.4 <- Reduce(rbind,My.Results.3) %>% na.omit()
  
  My.Results.4$Score.order <- My.Results.4$Score.order %>% as.integer() %>% as.character()
  
  My.Results.4$Group <- My.Results.4$Group %>% as.character()
  
  My.Results.4$Color <- My.Results.4$Group %>% factor() %>% sapply(.,function(x)RColorBrewer::brewer.pal(4,"Set1")[x])
  
  My.Results.4$Hamming.Distance <- My.Results.4$Hamming.Distance %>% as.character()
  
  
  My.Results.4$Score.order.plot <- sapply(1:nrow(My.Results.4),function(i){
    
    r <- My.Results.4[i,]
    if(r$Group=="G1")return(as.numeric(r$Score.order)-0.30)
    if(r$Group=="G2")return(as.numeric(r$Score.order)-0.10)
    if(r$Group=="G9")return(as.numeric(r$Score.order)+0.10)
    if(r$Group=="G10")return(as.numeric(r$Score.order)+0.30)
    
  })
  
  
  My.Results.6 <- My.Results.4
  
  
  
  # R    B Group
  # 1 0.00 0.00     1
  # 2 0.05 0.00     2
  # 3 0.10 0.00     3
  # 4 0.15 0.00     4
  # 5 0.20 0.00     5
  # 6 0.00 0.05     6
  # 7 0.00 0.10     7
  # 8 0.00 0.15     8
  # 9 0.00 0.20     9
  
  #My.Results.6 <<- My.Results.6
  
  Data.2.1 <- 
    My.Results.6 %>% 
    dplyr::group_by(.,Group,Score.order.plot,Group,Color) %>%
    dplyr::summarise(Score.Mean=mean(as.numeric(Score)),Score.SE=std(as.numeric(Score)))
  
  
  # plot 4 figures 
  
  # score order 
  p.2.1 <-
    ggplot(Data.2.1)+
    geom_errorbar(aes(x=Score.order.plot,
                      ymax=Score.Mean+Score.SE,
                      ymin=Score.Mean-Score.SE,
                      color= I(Color)),
                  width=0.2,size=0.5)+
    geom_col(aes(x=Score.order.plot,
                 y=Score.Mean,
                 fill= I(Color)),
             width=0.2)+
    scale_x_continuous(labels=c(1:10),breaks = c(1:10))+
    xlab("Rank of DELTA score")+
    ylab("DELTA Score")+
    theme_classic()
  
  
  # p value
  
  Data.2.2 <-
    My.Results.6 %>%
    dplyr::filter(Group %in% c("G1","G2","G9","G10"))%>%
    dplyr::mutate(Pvalue=-log10(as.numeric(pvalue))) %>%
    dplyr::mutate(Pvalue=ifelse(Pvalue>300,300,Pvalue)) %>% 
    dplyr::group_by(Group,Score.order.plot,Color) %>%
    dplyr::summarise(Pvalue.Mean=mean(as.numeric(Pvalue)),Pvalue.SE=std(Pvalue))
  
  
  p.2.2 <-
    ggplot(Data.2.2) +
    geom_errorbar(
      aes(
        x = Score.order.plot,
        ymax = Pvalue.Mean + Pvalue.SE,
        ymin = Pvalue.Mean - Pvalue.SE,
        color = I(Color)
      ),
      width=0.2,size=0.5,
      show.legend = F
    ) +
    geom_col(
      aes(x = Score.order.plot,
          y = Pvalue.Mean,
          fill = I(Color)),
      width = 0.2,
      show.legend = F
    )+
    scale_x_continuous(labels=c(1:10),breaks = c(1:10))+
    xlab("Rank of DELTA score")+
    ylab(expression(-log10(italic(P))))+
    theme_classic()
  #+geom_text(aes(x=Score.order.plot,
  #               y=Pvalue.Mean,
  #               label=Pvalue.Mean,
  #               color=I(Color)),size=1,show.legend = F)+
  
  
  # Tip size
  Data.2.3 <-
    My.Results.6 %>% 
    dplyr::filter(Group %in% c("G1","G2","G9","G10"))%>% 
    dplyr::group_by(Group,Score.order.plot,Color) %>% 
    dplyr::summarise(TipSize.Mean=mean(as.numeric(TipSize)),TipSize.SE=std(TipSize))
  
  
  p.2.3<-
    ggplot(Data.2.3)+
    geom_errorbar(aes(x=Score.order.plot,
                      ymax=TipSize.Mean+TipSize.SE,
                      ymin=TipSize.Mean-TipSize.SE,
                      color=I(Color)),
                  width=0.2,size=0.5
    )+
    geom_col(aes(x=Score.order.plot,
                 y=TipSize.Mean,
                 fill=I(Color)),
             width= 0.2
    )+
    scale_x_continuous(labels=c(1:10),breaks = c(1:10))+
    xlab("Rank of DELTA score")+
    ylab("Tip size")+
    theme_classic()
  
  
  
  # Hamming distance
  Data.2.4 <-
    My.Results.6 %>% 
    dplyr::filter(Group %in% c("G1","G2","G9","G10"))%>% 
    dplyr::group_by(Group,Score.order.plot,Color) %>% 
    dplyr::summarise(HD.Mean=mean(as.numeric(Hamming.Distance)),HD.SE=std(as.numeric(Hamming.Distance)))
  
  
  
  p.2.4 <-
    ggplot(Data.2.4)+
    geom_errorbar(aes(x=Score.order.plot,
                      ymax=HD.Mean+HD.SE,
                      ymin=HD.Mean-HD.SE,
                      color= I(Color),group=Group),
                  width=0.2,size=0.5
    )+
    geom_col(aes(x=Score.order.plot,
                 y=HD.Mean,
                 fill=I(Color),group=Group),
             width = 0.2)+
    scale_x_continuous(labels=c(1:10),breaks = c(1:10))+
    xlab("Rank of DELTA score")+
    ylab("Difference in expression status")+
    theme_classic()
  
  # legend
  p.legend <-
    get_legend(p.2.4 +
                 scale_color_manual(
                   labels = c(
                     "Full autonomous tree",
                     paste0(100*D1,"% terminal cells loss"),
                     "5% non-autonomous expression",
                     paste0(100*D1,"% terminal cell loss + \n5% non-autonomous expression")
                   ),
                   name = "Tree simulation model:",
                   values = I(RColorBrewer::brewer.pal(4, "Set1")[c(1, 3, 4, 2)])
                 ))
  
  
  p2 <-  cowplot::plot_grid(plotlist = list(p.2.1,p.2.3,p.2.2,p.2.4,p.legend),nrow=1)#,labels = c("C","D","E","F"))
  p2

}



# Calculate under each conditions
for(D1 in c(0.05,0.10,0.20,0.50,0.90)){
  
  score="common"
  m=1
  u=1
  prune=20
  test.num=100
  method="l"
  setwd("/mnt/data/home/yuanmeng/2017-2018/tree/SimulateTree/")
  dir.create(paste0("./data.common.","m",m,"u",u,"pr",prune,".R",D1,"B0.05"))
  setwd(paste0("./data.common.","m",m,"u",u,"pr",prune,".R",D1,"B0.05"))
  p2 <- MyPlot(score = "common",m=m,u=u,prune = prune,Cal = F,test.num=test.num,D1=D1)
  dir.create(paste0("../result2/",D1))
  figure_file <- paste0("../result2/",D1,"/Fig2B-E-",score,"-m",m,"u",u,"pr",prune,".pdf")
  ggsave(figure_file,p2,width = 25,height=3,units = "in")
  cat('Done!\n')
  setwd("../")
  
}











