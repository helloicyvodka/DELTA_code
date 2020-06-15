# outfile.folder <- "./test/"

F.ReadResult <- function(outfile.folder){
  
  outfile.list <- list.files(outfile.folder)
  
  

  r.list <-
    
    base::lapply(outfile.list,function(outfile){
      
    if(length(list.files(paste0(outfile.folder,outfile)))==0)return(NA)
    
    file <- paste0(outfile.folder,"/",outfile,"/",outfile,".alml")
    
    file.AllLin.df <- paste0(outfile.folder,"/",outfile,"/",outfile,".ReLin.All")
    
    file.cell.tree.all <- paste0(outfile.folder,"/",outfile,"/",outfile,".cell.tree.all")
    
    #file.LeafLin.df <-  paste0(outfile.folder,"/",outfile,"/",outfile,".ReLin.Leaf")
    
    
    
    AllLin.df <- read.table(file.AllLin.df,header = T,colClasses = "character",stringsAsFactors = F)
    
    #LeafLin.df <- read.table(file.LeafLin.df,header = T,colClasses = "character",stringsAsFactors = F)
    

    
    cell.tree.all.df <- read.table(file.cell.tree.all , header = T, colClasses = "character",stringsAsFactors = F)
    
    
    r <- ggvita::readal.alml(file)
    
    
    r.df <- data.frame(Info= rep(outfile,length(r)),Score.order=1:length(r),stringsAsFactors = F)
    
    r.df$Score <- list.mapv(r,.$Score)
    
    r.df$RootS <- list.mapv(r,.$RootS)
    
    r.df$RootT <- list.mapv(r,.$RootT)
    
    r.df$TipSize <- list.mapv(r,c(.$RootS,.$MatchS) %>% Find.tips() %>% length())
    
    r.df$pvalue <- list.mapv(r,.$PValue$pvalue)
    
    Info.df <- r.df$Info %>% lapply(function(x){ 
      t <-strsplit(x,split = "[a-zA-Z]") %>% unlist()
      t <- t[c(-1,-2)]
      t(as.matrix(t))
    }) %>% Reduce(rbind,.) %>% data.frame(stringsAsFactors = F) 
    
    colnames(Info.df ) <- c("ID","a","N","K","tmax","Dmax","bias","RdropRate")
    
    r.df <- cbind(r.df,Info.df)
    
    r.df <- data.table(as.matrix(r.df))
    
    r.df <- data.frame(r.df)
    
    r.df$Gene.OF.RootS <- sapply(1:nrow(r.df),function(i){
      
      rr <- r.df[i,]
      rr.RootS <- rr["RootS"] %>% as.character()
      
      #rr.RootT <- rr["RootT"] %>% as.character()
      
      rr.RootS.Ori <- AllLin.df[AllLin.df$NewAllLin==rr.RootS,"AllLin"] %>% as.character()
      #rr.RootT.Ori <- AllLin.df[AllLin.df$NewAllLin==rr.RootT,"AllLin"] %>% as.character()
      rr.RootS.Ori <- ifelse(rr.RootS.Ori=="Root","",rr.RootS.Ori)
      
      Gene.OF.RootS <- cell.tree.all.df[cell.tree.all.df$Lineage==rr.RootS.Ori,"Gene.OF"] %>% as.character()
      #Gene.OF.RootT <- cell.tree.all.df[cell.tree.all.df$Lineage==rr.RootT.Ori,"Gene.OF"] %>% as.character()
      
      Gene.OF.RootS 
      
    })
    
    r.df$Gene.OF.RootT <- sapply(1:nrow(r.df),function(i){
      
      
      rr <- r.df[i,]
      #rr.RootS <- rr["RootS"] %>% as.character()
      rr.RootT <- rr["RootT"] %>% as.character()
      
      #rr.RootS.Ori <- AllLin.df[AllLin.df$NewAllLin==rr.RootS,"AllLin"] %>% as.character()
      rr.RootT.Ori <- AllLin.df[AllLin.df$NewAllLin==rr.RootT,"AllLin"] %>% as.character()
      
      rr.RootT.Ori<- ifelse(rr.RootT.Ori=="Root","",rr.RootT.Ori)
      
      #Gene.OF.RootS <- cell.tree.all.df[cell.tree.all.df$Lineage==rr.RootS.Ori,"Gene.OF"] %>% as.character()
      Gene.OF.RootT <- cell.tree.all.df[cell.tree.all.df$Lineage==rr.RootT.Ori,"Gene.OF"] %>% as.character()
      
      Gene.OF.RootT 
      
    })
    
    
    r.df$Hamming.Distance <- apply(r.df,1,function(x) sum(unlist(strsplit(x["Gene.OF.RootS"],"")) != unlist(strsplit(x["Gene.OF.RootT"],""))))
    
    
    
    
    r.df
    
  }) 
  
  r.list %>% Reduce(rbind,.) 

  
}




  
  



