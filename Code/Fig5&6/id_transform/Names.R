x<-readLines("./id_transform/worm_proteinID_genename.txt")


x<-x[lapply(x,function(s){
  p<-strsplit(s,split = " ") %>% unlist() %>% length()
  p==2
  }) %>% unlist()]
writeLines(x,"./id_transform/worm_proteinID_genename2.txt")