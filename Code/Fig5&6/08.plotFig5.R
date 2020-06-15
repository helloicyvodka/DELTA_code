

require(dplyr)
require(ggplot2)
require(ggtree)
require(ggvita)

hom.trans.Score_order %>% filter(Class.digital=="Yes") %>% arrange(-Score) %>% View



My.Fig <-
function(gene,score.order){
  
  #tmp.data <-hom.trans.result.common[[paste0(gene,"_.almg")]]$mut_wt[[as.character(score.order)]]
  tmp.p <-
    ggvita::ggvita(hom.trans.result.common[[paste0(gene,"_.almg")]]$mut_wt,score.order+1) %++%
    ggtree::geom_tippoint(ggplot2::aes(fill = I(tip.fill)), size = 2, shape = 21, color = "NA")
  
  tmp.p%++% ggvita::stat_prune(tmp.p)
  
}

My.Fig.5 <- function(gene,score.order,...){
  
  ttt <-
  ggvita::ggvita(All.gene.alml.list.2[[gene]],score.order,...) %++% 
    ggtree::geom_tippoint(ggplot2::aes(fill = I(tip.fill)), size = 2, shape = 21, color = "NA")
  
  ttt %++% stat_prune(ttt)
  
}



# digital detected

My.Fig.5("MOM-2",6)
My.Fig.5("CUL-1",12)



tmp_pp <- My.Fig.5("CPSF-1",29)
tmp_pp %++% geom_tiplab(aes(label=node.seq),size=1)


pdf("./class_digital.pdf",paper = "a4",onefile = T)
df <- hom.trans.Score_order %>% filter(Class.digital=="Yes")
ttt <- lapply(1:nrow(df),function(i){
  
  r <- df[i,]
  My.Fig(as.character(r$Gene),r$Score_order.mut_wt+1)%++%
    geom_text(x=3,y=5,label=i)
  
})
ttt
dev.off()



# 15 22 26 27 28 35 42



pdf("./class_DELTA2.pdf",paper = "a4",onefile = T)
df <- hom.trans.Score_order %>% dplyr::filter(Class.digital=="No",Score_order.mut_wt<20,Score.order.elevation>20)
ttt <- parallel::mclapply(1:100,function(i){
  
  r <- df[i,]
  My.Fig(as.character(r$Gene),r$Score_order.mut_wt+1)%++%
    ggplot2::geom_text(x=3,y=5,label=i)
  
},mc.cores=20)
ttt
dev.off()

# 24 28 31 42 60 71
# 28


i=28



My.Fig("MOM-2",6)
My.Fig("CAMT-1",4)



