
# Hyper parameter -------------------------------------

D1=0.05
parent.folder="./"
B <- c(0.05,0.05,0.05,0.05)
D <- c(D1,0.40,0.45,0.50)
Parm.df <- data.frame(R=c(D,0,0,0,0),B=c(0,0,0,0,B),stringsAsFactors = F)
Parm.df <- rbind(c(0,0),Parm.df,c(D1,0.05))
Parm.df$Group <- sapply(as.character(1:nrow(Parm.df)),function(x)paste0("G",x))

#Get the correlation heatmap --------------------------


Get.cor.df <- function(R,B,parent.folder="./"){
  
  outputfolder <- paste0(parent.folder,"data.R",R,"B",B,"/")    #./data.R0B0/
  
  All.Try.folder <- list.files(outputfolder,full.names = T)
  All.Try.folder.Info <- list.files(outputfolder,full.names = F)
  
  
  lapply(All.Try.folder.Info,function(f){
    
    outputfolder.one <- paste0(outputfolder,f)  
    one.alml <- paste0(outputfolder.one,"/",f,".alml")
    one.alml.result <- ggvita::readal.alml(one.alml)
    one.cell.tree.all.1 <- paste0(outputfolder.one,"/",f,".1.cell.tree.all")
    one.cell.tree.all.2 <- paste0(outputfolder.one,"/",f,".2.cell.tree.all")
    
    one.cell.tree.all.1.df <- 
      read.table(one.cell.tree.all.1,header = T, colClasses = "character",stringsAsFactors = F)
    one.cell.tree.all.2.df <- 
      read.table(one.cell.tree.all.2,header = T, colClasses = "character",stringsAsFactors = F)
    
    one.cell.tree.all.1.df$Lineage[one.cell.tree.all.1.df$Lineage==""] <- "Root"
    one.cell.tree.all.2.df$Lineage[one.cell.tree.all.2.df$Lineage==""] <- "Root"
    
    Expr.list <-
      lapply(1:10,function(s.o){
        
        Lin.S <- c(one.alml.result[[s.o]]$RootS,one.alml.result[[s.o]]$MatchS)
        Lin.T <- c(one.alml.result[[s.o]]$RootT,one.alml.result[[s.o]]$MatchT)
        if(R!=0){
          
          S.ReLin.path <- paste0(outputfolder.one,"/",f,".1.ReLin.All")
          T.ReLin.path <- paste0(outputfolder.one,"/",f,".2.ReLin.All")
          
          S.ReLin.df <-     
            read.table(S.ReLin.path,header = T, colClasses = "character",stringsAsFactors = F)
          T.ReLin.df <-     
            read.table(T.ReLin.path,header = T, colClasses = "character",stringsAsFactors = F)
          
          S.ReLin <- S.ReLin.df$AllLin
          names(S.ReLin) <- S.ReLin.df$NewAllLin
          
          T.ReLin <- T.ReLin.df$AllLin
          names(T.ReLin) <- T.ReLin.df$NewAllLin
          
          Lin.S <- S.ReLin[Lin.S] %>% as.character()
          Lin.T <- T.ReLin[Lin.T] %>% as.character()
          
          
          
        }
        
        Lin.S.expr <- sapply(Lin.S,function(cell){
          
          gene.of <- one.cell.tree.all.1.df[one.cell.tree.all.1.df$Lineage==cell,"Gene.OF"]
          ifelse(is.null(gene.of),NA,gene.of)
          
        },USE.NAMES = F)
        
        Lin.T.expr <- sapply(Lin.T,function(cell){
          
          gene.of <- one.cell.tree.all.2.df[one.cell.tree.all.2.df$Lineage==cell,"Gene.OF"]
          ifelse(is.null(gene.of),NA,gene.of)
          
        },USE.NAMES = F)
        
        Expr.df <- data.frame(Score.order=s.o,Lin.S=Lin.S,Lin.T=Lin.T,Lin.S.expr=Lin.S.expr,Lin.T.expr=Lin.T.expr,stringsAsFactors = F)
        Expr.df$Lin.S.expr <- Expr.df$Lin.S.expr %>% strsplit(split = "") %>% lapply(.,as.numeric,USE.NAMES = F)
        Expr.df$Lin.T.expr <- Expr.df$Lin.T.expr %>% strsplit(split = "") %>% lapply(.,as.numeric,USE.NAMES = F)
        
        Expr.cor <-  sapply(1:N,function(i){
          e.S <- sapply(Expr.df$Lin.S.expr,function(e)e[i])
          e.T <- sapply(Expr.df$Lin.T.expr,function(e)e[i])
          cor(e.S,e.T,method = "pearson")
        },USE.NAMES = F)
        
        Expr.cor
      })  
    
    Expr.ld <- Reduce(cbind,Expr.list)
    colnames(Expr.ld) <- as.character(1:10)
    #sapply(1:10,function(x)paste0("Order",x))
    rownames(Expr.ld) <- as.character(1:8)
    #sapply(1:8,function(x)paste0("Gene",x))
    Expr.ld2 <- reshape2::melt(Expr.ld)
    colnames(Expr.ld2) <- c("Gene","Score.order","Cor")
    Expr.ld2
    
  })
  
}








Cor.plot.list <- list()
for(i in c(9) ){
  r <- Parm.df[i,]
  df.list <- Get.cor.df(r$R,r$B)
  df <- Reduce(rbind,df.list)
  df2 <- df %>% group_by(Gene,Score.order) %>% summarise(Cor.Mean=mean(Cor,na.rm = T))
  Cor.plot.list[[i]] <- df2
  cat(i,"\n")
}

for(i in c(1,2,9,10)){
  
  Cor.plot.list[[i]]$Group <- i
  
}


Cor.plot.list.df <- Reduce(rbind,Cor.plot.list)


ggplot(Cor.plot.list.df,aes(group=Group))+
  geom_tile(aes(x=factor(Score.order),y=factor(Gene),fill=Cor.Mean),color="white",size=0.2)+
  theme_gray()+
  xlab("Score order")+
  ylab("Gene")+
  theme(axis.text = element_text(size=12),axis.ticks = element_blank(),panel.background = element_blank())+
  scale_fill_gradientn(colors = rev(rainbow(7)))+
  facet_wrap(c("Group"),ncol=2)






D1=0.90
parent.folder="./"


B <- c(0.05,0.05,0.05,0.05)

D <- c(D1,0.40,0.45,0.50)



Parm.df <- data.frame(R=c(D,0,0,0,0),B=c(0,0,0,0,B),stringsAsFactors = F)

Parm.df <- rbind(c(0,0),Parm.df,c(D1,0.05))

Parm.df$Group <- sapply(as.character(1:nrow(Parm.df)),function(x)paste0("G",x))






# Get hamming distance heatmap ---------------------------------------------------------


Get.hd.df <- function(R,B,parent.folder="./"){
  
  outputfolder <- paste0(parent.folder,"data.R",R,"B",B,"/")    #./data.R0B0/
  
  All.Try.folder <- list.files(outputfolder,full.names = T)
  All.Try.folder.Info <- list.files(outputfolder,full.names = F)
  
  
  lapply(All.Try.folder.Info,function(f){
    
    outputfolder.one <- paste0(outputfolder,f)  
    one.alml <- paste0(outputfolder.one,"/",f,".alml")
    if(file.exists(one.alml)){
      one.alml.result <- ggvita::readal.alml(one.alml)  
    }else{return(NA)}
    
    one.cell.tree.all.1 <- paste0(outputfolder.one,"/",f,".1.cell.tree.all")
    one.cell.tree.all.2 <- paste0(outputfolder.one,"/",f,".2.cell.tree.all")
    
    one.cell.tree.all.1.df <- 
      read.table(one.cell.tree.all.1,header = T, colClasses = "character",stringsAsFactors = F)
    one.cell.tree.all.2.df <- 
      read.table(one.cell.tree.all.2,header = T, colClasses = "character",stringsAsFactors = F)
    
    one.cell.tree.all.1.df$Lineage[one.cell.tree.all.1.df$Lineage==""] <- "Root"
    one.cell.tree.all.2.df$Lineage[one.cell.tree.all.2.df$Lineage==""] <- "Root"
    
    Expr.list <-
      lapply(1:10,function(s.o){
        
        Lin.S <- c(one.alml.result[[s.o]]$RootS,one.alml.result[[s.o]]$MatchS)
        Lin.T <- c(one.alml.result[[s.o]]$RootT,one.alml.result[[s.o]]$MatchT)
        if(R!=0){
          
          S.ReLin.path <- paste0(outputfolder.one,"/",f,".1.ReLin.All")
          T.ReLin.path <- paste0(outputfolder.one,"/",f,".2.ReLin.All")
          
          S.ReLin.df <-     
            read.table(S.ReLin.path,header = T, colClasses = "character",stringsAsFactors = F)
          T.ReLin.df <-     
            read.table(T.ReLin.path,header = T, colClasses = "character",stringsAsFactors = F)
          
          S.ReLin <- S.ReLin.df$AllLin
          names(S.ReLin) <- S.ReLin.df$NewAllLin
          
          T.ReLin <- T.ReLin.df$AllLin
          names(T.ReLin) <- T.ReLin.df$NewAllLin
          
          Lin.S <- S.ReLin[Lin.S] %>% as.character()
          Lin.T <- T.ReLin[Lin.T] %>% as.character()
          
          
          
        }
        
        Lin.S.expr <- sapply(Lin.S,function(cell){
          
          gene.of <- one.cell.tree.all.1.df[one.cell.tree.all.1.df$Lineage==cell,"Gene.OF"]
          ifelse(is.null(gene.of),NA,gene.of)
          
        },USE.NAMES = F)
        
        Lin.T.expr <- sapply(Lin.T,function(cell){
          
          gene.of <- one.cell.tree.all.2.df[one.cell.tree.all.2.df$Lineage==cell,"Gene.OF"]
          ifelse(is.null(gene.of),NA,gene.of)
          
        },USE.NAMES = F)
        
        Expr.df <- data.frame(Score.order=s.o,Lin.S=Lin.S,Lin.T=Lin.T,Lin.S.expr=Lin.S.expr,Lin.T.expr=Lin.T.expr,stringsAsFactors = F)
        Expr.df$Lin.S.expr <- Expr.df$Lin.S.expr %>% strsplit(split = "") %>% lapply(.,as.numeric,USE.NAMES = F)
        Expr.df$Lin.T.expr <- Expr.df$Lin.T.expr %>% strsplit(split = "") %>% lapply(.,as.numeric,USE.NAMES = F)
        
        #Expr.df$Hamming.distance <- apply(Expr.df,1,function(x)sum(x$Lin.S.expr!= x$Lin.T.expr))
        
        N=8
        
        Expr.HD <-  sapply(1:N,function(i){
          e.S <- sapply(Expr.df$Lin.S.expr,function(e)e[i])
          e.T <- sapply(Expr.df$Lin.T.expr,function(e)e[i])
          sum(e.S!=e.T)/length(e.S)
        },USE.NAMES = F)
        
        Expr.HD
      })  
    
    Expr.ld <- Reduce(cbind,Expr.list)
    colnames(Expr.ld) <- as.character(1:10)
    #sapply(1:10,function(x)paste0("Order",x))
    rownames(Expr.ld) <- as.character(1:8)
    #sapply(1:8,function(x)paste0("Gene",x))
    Expr.ld2 <- reshape2::melt(Expr.ld)
    colnames(Expr.ld2) <- c("Gene","Score.order","HD")
    Expr.ld2
    
  })
  
}


HD.plot.list <- list()
for(i in c(1,2,9,10) ){
  r <- Parm.df[i,]
  df.list <- Get.hd.df(r$R,r$B)
  df <- Reduce(rbind,df.list)
  df2 <- df %>% dplyr::group_by(Gene,Score.order) %>%dplyr::summarise(HD.Mean=mean(HD,na.rm = T))
  HD.plot.list[[i]] <- df2
  cat(i,"\n")
}


for(i in c(1,2,9,10)){
  HD.plot.list[[i]]$Group <- i
}


HD.plot.list.df <- Reduce(rbind,HD.plot.list)

HD.plot.list.df$Group <- factor(HD.plot.list.df$Group)
levels(HD.plot.list.df$Group) <- c("Full autonomous tree",
                                   paste0(D1*100,"% terminal cells loss"),
                                   "5% non-autonomous expression",
                                   paste0(D1*100,"% terminal cell loss + \n5% non-autonomous expression"))


ggplot(HD.plot.list.df %>% filter(! Gene %in% c(1)) %>% filter(!is.na(Gene)),aes(group=Group))+
  geom_tile(aes(x=factor(Score.order),y=factor(Gene),fill=HD.Mean),color="white",size=0.2)+
  theme_gray()+
  xlab("Score order")+
  ylab("Gene")+
  theme(axis.text = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()
  )+
  scale_fill_gradient2(low = "navy",
                       mid = "gray",
                       high= "red",
                         midpoint = 0.25,
                         #colors = rev(RColorBrewer::brewer.pal(7,"RdYlBu")),
                         name="Normalized Hamming distance",
                       limits=c(0,0.5)
                       )+
  facet_wrap(c("Group"),ncol=2)




