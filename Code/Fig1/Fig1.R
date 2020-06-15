
library(magrittr)
library(ggvita)
library(ggplot2)
library(dplyr)
setwd("~/2017-2018/tree/C.elegans_VS_C.elegans")



# create and import cost files --------------------------------------------
cost.origin <- read.table("./cost.tsv",header = F,stringsAsFactors = F)
cost.YXJ<- read.table("./matrix_fun.csv",header = F,stringsAsFactors = F)
cost.m8u8<- data.frame(cost.origin[,1:2],stringsAsFactors = F)
cost.m8u8 <- cost.m8u8 %>% mutate(V3=ifelse(V1==V2,8,-8))
write.table(cost.m8u8 ,"./cost.common.m8u8.tsv",row.names = F,col.names = F,quote = F)


# hyper parameters --------------------------------------------------------
cost <- c("./matrix_fun.csv",
          "./cost.common.m8u8.tsv",
          "./cost.tsv",
          "/mnt/data/home/ATPs/W/Celegans/Lineages/costMatch10otherMinus100.tsv")[1]
outfile <- "./CXL-pr1-3000.alml"

fileS="~/2017-2018/tree/C.elegans_VS_C.elegans/data/fun.alm"
fileT = fileS
max_target = 10000
pr=100
test = 10


# Calculating cell types matching -----------------------------------------

C.vs.C <- ggvita::readal(outfile = outfile,
                         method = "g",
                         fileS = fileS,
                         fileT = fileT,
                         cost = cost,
                         test = test,
                         prune = pr)


cell_dict <- c("Str"="Structural","Neu"="Neuron","Mus"="Muscle","Int"="Intestine","Gla"="Gland","Ger"="Germ","Epi"="Epithelial","Dea"="Death","Bla"="Blast")
CC_p <- ggvita(C.vs.C)
CC_p_2 <- CC_p$plot$ggS$data$Class %>% table() %>% data.frame()
colnames(CC_p_2) <- c("CellTypes","Count")
CC_p_2$CellTypes <- sapply(CC_p_2$CellTypes,function(x)cell_dict[x])
CC_p_2$CellTypes2 <- CC_p_2$CellTypes

# Plot cell types matching -----------------------------------------
ggplot(CC_p_2)+
  geom_tile(aes(x=`CellTypes`,y=`CellTypes2`,fill=Count))+
  geom_text(aes(x=`CellTypes`,y=`CellTypes2`,label=Count))+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9,"YlOrRd")[c(2,3,4,5,7,8,9)])+
  ylab("Isomorphic C. elegans CLT")+
  xlab("C. elegans CLT")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=-30,hjust=0.02),
        axis.ticks = element_blank())



