



# plot multiple figures to find homeotic transformation examples

data2 <-
  hom.trans.Score_order %>% 
  filter(is.na(Elevation.order)!=T) 

data3 <- 
  hom.trans.Score_order %>% 
  filter(Class.digital=="Yes")

p <- 
  ggplot(fortify(data2)) + 
  geom_tile(aes(x=as.integer(Elevation.order),y=Gene,fill = log10(as.numeric(Score.mut_wt))),color =NA) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(8,"Greys"),name="log10(Score)")+
  labs(title="Homoetic transformation")+
  theme(axis.text.y = element_text(size = 3),plot.title = element_text(hjust = 0.5))+
  ylab("Gene")+
  xlab("Order of score order elevation ")+
  scale_x_continuous(breaks = c(0:10)*100)+
  geom_tile(data=fortify(data3),aes(as.integer(Elevation.order),Gene),color="red",fill="red")

p

ggsave("./result2/100_Find.pdf",width = 8,height = 8,units = "in")



# plot a single example to see

My.Fig <- function(gene,score.order){
  
  #tmp.data <-hom.trans.result.common[[paste0(gene,"_.almg")]]$mut_wt[[as.character(score.order)]]
  tmp.p <-
    ggvita::ggvita(hom.trans.result[[gene]]$mut_wt,score.order) %++%
    ggtree::geom_tippoint(ggplot2::aes(fill = I(tip.fill)), size = 2, shape = 21, color = "NA")
  
  tmp.p%++% ggvita::stat_prune(tmp.p)
  
}


My.Fig("ABI-1",11)

hom.trans.Score_order %>% 
  dplyr::filter(RootS %in% Fate_change.s3$Alternative.Path.from 
                & RootT %in% Fate_change.s3$Alternative.Path.to) %>% 
  dplyr::
  
  
  hom.trans.Score_order.2 <- hom.trans.Score_order
hom.trans.Score_order.2$RootS[hom.trans.Score_order.2$RootS=="Root"] <- ""
hom.trans.Score_order.2$RootT[hom.trans.Score_order.2$RootT=="Root"] <- ""
hom.trans.Score_order.2$Nchar.mean<- mclapply(1:nrow(hom.trans.Score_order.2),function(x){
  
  r <- hom.trans.Score_order.2[x,]
  
  mean(nchar(r$RootS),nchar(r$RootT))
  
  
  
},mc.cores = 20) %>% unlist()


ggplot(hom.trans.Score_order.2)+
  stat_bin(aes(x=Score_order,fill=Class.sister,alpha=-log(Nchar.mean)),binwidth =10)+
  labs(title="Symmeteric pairs")+
  xlab("Score order")+
  theme(plot.title = element_text(hjust =0.5))+
  scale_fill_manual(values=c("skyblue","red"),name="Symmeteric pairs")


head(hom.trans.Score_order.2 %>% 
       dplyr::filter(Class.digital=="Yes") %>% 
       dplyr::arrange(-as.numeric(Score)),
     n=10
)




# 
#       Gene RootS RootT Score Score_order Class.digital Class.sister Nchar.mean
# 1    WWP-1    00    01  3521           2           Yes           No          2
# 2    APX-1    01    00  3084           2           Yes           No          2
# 3    LAG-1    01    00  2632           4           Yes           No          2
# 4    MOM-2   001   000  1854           6           Yes           No          3
# 5    DIV-1     0    01  1702           3           Yes           No          1
# 6    NPP-4    00    01  1523           2           Yes           No          2
# 7    SEL-8    01    00  1383           5           Yes           No          2
# 8  F23F1.5    00    01  1286           4           Yes           No          2
# 9    PAR-6    11    10  1159          16           Yes           No          2
# 10  CEH-24    00    01  1146           4           Yes           No          2


outfile <- "./data3/mut_wt_local/WWP-1_.alml"
fileS <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/mut-tree-tips/WWP-1_.alm"   
fileT <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/wt-tree-tips/WWP-1_.alm"   
cost <-  "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"


All.gene.alml.list <- 
  
  parallel::mclapply(result.list.files.genes,function(x){
    
    gene <- as.character(x)
    
    outfile <-  paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/mut_wt_local/",gene,"_.alml")
    fileS <- paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/mut-tree-tips/",gene,"_.alm")
    fileT <- paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/wt-tree-tips/",gene,"_.alm")
    cost <-  "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"
    
    readal(outfile,fileS,fileT,cost,method = "l")
    
  },mc.cores = 20)



names(All.gene.alml.list) <- result.list.files.genes





All.gene.alml.list.2 <- 
  
  parallel::mclapply(result.list.files.genes,function(x){
    
    gene <- as.character(x)
    
    outfile <-  paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data3/common-alml/mut-vs-wt/",gene,"_.alml")
    fileS <- paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/mut-tree-tips/",gene,"_.alm")
    fileT <- paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data2/wt-tree-tips/",gene,"_.alm")
    cost <-  "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"
    
    readal(outfile,fileS,fileT,cost,method = "l")
    
  },mc.cores = 20)



names(All.gene.alml.list.2) <- result.list.files.genes