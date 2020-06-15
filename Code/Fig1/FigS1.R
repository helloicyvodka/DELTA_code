library(dplyr)
library(purrr)
library(ggplot2)
library(ggtree)

# create isomophic celegans tree

fun.alm <- read.table("./fun.alm",header = T,stringsAsFactors = F,colClasses = "character")

fun.alm2 <- rbind(fun.alm,data.frame(Lineage=fun.alm$Lineage %>% Find.tree.from.tips() %>% setdiff(fun.alm$Lineage),
                                     Name=NA,Class=NA))

fun.alm2 <- fun.alm2 %>% arrange(Lineage)

fun.alm2 <- fun.alm2[fun.alm2$Lineage!="",]

fun.alm2$Lineage.Sister <-
  fun.alm2$Lineage %>% sapply(function(x){
  paste0(substr(x,1,nchar(x)-1),
         ifelse(substr(x,nchar(x),nchar(x))=="0","1","0"))

}) %>% as.character()

fun.alm2$If.Filp <- runif(nrow(fun.alm2),0,1)<0.3

fun.alm3 <- fun.alm2 %>% filter(If.Filp==T)

for(i in 1:nrow(fun.alm3)){
  cat(i,"\n")
  r <- fun.alm3[i,]

  fun.alm2$Lineage<- fun.alm2$Lineage %>% sapply(function(ii){

    if(startsWith(ii,prefix = r$Lineage))sub(pattern = r$Lineage,replacement = r$Lineage.Sister,ii)
    else if(startsWith(ii,prefix = r$Lineage.Sister))sub(pattern = r$Lineage.Sister,replacement = r$Lineage,ii)
    else(ii)
  }) %>% as.character()

}

write.alm(fun.alm2[fun.alm2$Lineage %in% Find.tips(c("",fun.alm2$Lineage)),c(1,2,3)],file.name = "~/2017-2018/tree/ReDrawResult/fun2.alm")



# calculate isomophic celegans and original celegans CLTs
DELTA(treeS = "~/2017-2018/tree/ReDrawResult/fun.alm",
      treeT = "~/2017-2018/tree/ReDrawResult/fun2.alm",
      cost = "~/2017-2018/tree/ReDrawResult/cost2.tsv",
      outfile = "~/2017-2018/tree/ReDrawResult/fun2.alml",
      method =  "l",
      max_target = 100,all = "T",prune = "1",test = "100",DELTA.address = DELTA.address)

C.ele.iso <- ggvita::readal(fileS = "~/2017-2018/tree/ReDrawResult/fun.alm",
               fileT = "~/2017-2018/tree/ReDrawResult/fun2.alm",
               cost = "~/2017-2018/tree/ReDrawResult/cost2.tsv",
               outfile = "~/2017-2018/tree/ReDrawResult/fun2.alml",
               method =  "l",
               all = T,prune = "1",test = "10")


p <-
  ggvita(alml_list = C.ele.iso,result.order = 1,trace_down_for_pruned = T) %++%
  geom_tippoint(aes(fill=I(tip.fill)),size=1,shape=21,color="NA")




pdf(file = "./Fig1F.pdf",paper = "a4")
DrawScoreMatrix(C.ele.iso)
dev.off()















