
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggvita)
library(data.table)

setwd("~/2017-2018/tree/Pmarina_Celegans")

ggvita::DELTA(treeS = "./data/pma-new.alm",
              treeT = "./data/fun.alm",
              cost  =  "./data/matrix_vs_f_fun.csv",
              outfile ="./data/PvsC3.almg",
              method = "g",
              prune = 80,
              test =  100,
              all = F)




P.vs.C <- ggvita::readal(outfile = "./data/PvsC3.almg",
                         fileS = "./data/pma-new.alm",
                         fileT = "./data/fun.alm",
                         cost = "./data/matrix_vs_f_fun.csv",
                         method = "g",all = "F",prune = 80)




p <-
  ggvita(alml_list = P.vs.C,
         result.order = 1,
         size=0.1,
         trace_down_for_pruned = T) %++%
  geom_tippoint(aes(fill=Class),size=3,shape=21,color="NA",show.legend = T)
#p


p$plot$ggS$data$Class[p$plot$ggS$data$Class=="Epi"] <- "Epi.P"
p$plot$ggT$data$Class[p$plot$ggT$data$Class=="Epi"] <- "Epi.C"

p <-
  p %++%
  ggplot2::scale_fill_manual(values = c("X" = "black", "Dea" = "black",
                                        "Ner" = "aquamarine","Neu" = "aquamarine1", "Str" = "aquamarine3",
                                        "Pha" = "gold3", "Gla" = "gold2",
                                        "Mus" = "chocolate",
                                        "Epi.P" = "skyblue", "Epi.C" = "skyblue2",
                                        "?" = "grey", "Bla" = "grey60" ,
                                        "Int" = "hotpink",
                                        "Ger" = "purple"),name="Class")%++%
  ggtree::theme(legend.position = "right",legend.background = ggplot2::element_rect(size = 0.8))



#pdf(file = "../Pmarina VS Celegans/20181210-2.pdf",height = 5,width = 72)

pp <-
  p%++%
  ggvita::stat_prune(p,size =3,color = "blue")%++%
  ggtree::geom_tiplab(aes(label=node.seq.ori),size=0.5,hjust = 0.5)

#print(pp,show.legend = F)

#dev.off()




pma_new.alm <- data.table::fread("../Pmarina VS Celegans/data/pma-new.alm",colClasses = "character")
pma_new.alm$LN <- ggvita::LN_to_Bin(pma_new.alm$Name)
#pma_new.alm[!pma_new.alm$Lineage==pma_new.alm$LN,]



# Merge Pmarina and Celegans data ----------------------------------------



data.match <- data.table(MatchS=pp$data$MatchS,MatchT=pp$data$MatchT,stringsAsFactors = F)
data.S <- pp$plot$ggS$data
colnames(data.S) <- paste0("S_",colnames(data.S))
data.T <- pp$plot$ggT$data
colnames(data.T) <- paste0("T_",colnames(data.T))

data.match <- merge(data.match,data.S,by.x="MatchS",by.y="S_node.seq.ori",all.x = T)
data.match <- merge(data.match,data.T,by.x="MatchT",by.y="T_node.seq.ori",all.x = T)

data.match.cell_type <- data.match[S_isTip==T&T_isTip==T,.N,by=c("S_Class","T_Class")]

colnames(data.match.cell_type) <- c("P.marina","C.elegans","N")

cell_dict <- c("Epi.C"="Epi","X"="Dea","Epi.P"="Epi","?"="Other")

data.match.cell_type$P.marina <- sapply(data.match.cell_type$P.marina,function(i){
  if(i %in% names(cell_dict)) i <- cell_dict[i]
  return(i)
})
data.match.cell_type$C.elegans <- sapply(data.match.cell_type$C.elegans,function(i){
  if(i %in% names(cell_dict)) i <- cell_dict[i]
  return(i)
})


ggplot(data.match.cell_type)+
  geom_tile(aes(x=`P.marina`,y=`C.elegans`,fill=N))+
  geom_text(aes(x=`P.marina`,y=`C.elegans`,label=N))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9,"YlOrRd")[c(2,3,4,5,7,8,9)])+
  theme_bw()+
  theme(panel.grid = element_blank())



# Transform Lineage -------------------------------------------------------

source("~/2017-2018/tree/transform ppp to 111/Bin_to_LN.R")
#source("~/2017-2018/tree/transform ppp to 111/Bin_to_LN_pma.R")

data.match$MatchS.Lin <- data.match$MatchS %>% mclapply(Bin_to_LN,mc.cores=10) %>% unlist()
data.match$MatchT.Lin <- data.match$MatchT %>% mclapply(Bin_to_LN,mc.cores=10) %>% unlist()

#tmp <- data.match$MatchT %>% mclapply(Bin_to_LN,mc.cores=10)
#data.match$MatchT.Lin <- unlist(tmp)

tmp <- data.match$MatchS %>% sapply(Bin_to_LN_pma)

#data.match[nchar(MatchS.Lin)==0|nchar(MatchT.Lin)==0]
#data.match[length(MatchS.Lin)==0|length(MatchT.Lin)==0]

#data.match[MatchS.Lin==""|MatchS.Lin==""]


# for(i in 1:nrow(data.match)){
# 
#   r <- data.match[i,]
# 
#   if(length(r$MatchS.Lin[[1]])==0)data.match[i,"MatchS.Lin"] <- r$S_Name
#   if(length(r$MatchT.Lin[[1]])==0)data.match[i,"MatchT.Lin"] <- r$S_Name
# 
# }


fwrite(data.match[,c("MatchS.Lin","MatchT.Lin","S_isTip","T_isTip","S_Class","T_Class","S_Name","T_Name","S_mp","T_mp")],"./cell_type_match.csv",quote=T)


for(i in 1:length(data.match$MatchT)){

  print(data.match$MatchT[i])
  print(data.match$MatchT[i] %>% Bin_to_LN)
}


for(i in 1:length(data.match$MatchS)){
  
  print(data.match$MatchS[i])
  print(data.match$MatchS[i] %>% Bin_to_LN_pma)
}










