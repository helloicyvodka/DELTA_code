library(rlist)
library(parallel)
library(magrittr)
library(pipeR)
library(rlist)

# wt <- as.data.frame(fread("~/2017-2018/tree/tree_fate_change/data/data_from_digital/WT_merge/HYP_NEU_PHA.txt"))

# source("~/2017-2018/tree/ggvita/R/LN_to_Bin.R")

# wt$Lineage <- LN_to_Bin(wt$cell_name)

wt <- com.expr.2

wt.2 <- data.frame(Lineage=wt$Lineage,Name=wt$cell_name,Class=wt$Class.num)

write.table(wt.2,file="./data/data_from_digital/WT_merge/wt_full_for_DELTA_.txt",quote = F,row.names = F)

rm_tail_blank_line("../WT_merge/wt_full_for_DELTA.txt")



get_parent <- function(x) {
  x <- as.character(x)
  x <- substr(x, 1, (nchar(x) - 1))
  if (x == "") {
    x <- "Root"
  }
  x
}






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
  
  mut <- mut.full.df.list[[x]]
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








