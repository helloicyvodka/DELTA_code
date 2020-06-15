F.Simulate.Dropped.Leaf <- function(tree1,RdropRate){
  
  tree.dfLeaf <- ToDataFrameTable(tree1,"pathString","level","name",'Lineage', 'DAge','tAge',
                                  Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))
  
  tree.dfLeaf$Drop <- ifelse(runif(nrow(tree.dfLeaf),0,1)>=RdropRate,0,1)
  
  # tree.dfLeaf.drop.unpaired <- tree.dfLeaf %>% filter(Drop==1) %>% `$`("Lineage") %>% substr(1,nchar(.)-1) %>% table 
  # 
  # tree.dfLeaf.drop.unpaired <- names(tree.dfLeaf.drop.unpaired[tree.dfLeaf.drop.unpaired==1])
  # 
  # tree.dfLeaf$Lineage <- sapply(1:nrow(tree.dfLeaf),function(x){
  #   
  #   r <- tree.dfLeaf[x,]
  #   
  #   mother <- as.character(r["Lineage"]) %>% substr(1,nchar(.)-1)
  #   
  #   if( (mother %in% tree.dfLeaf.drop.unpaired ) &
  #       (r$Drop==0)
  #   ){
  #     
  #     mother
  #     
  #   }else{
  #     
  #     r$Lineage  
  #     
  #   }
  #   
  # })
  
  tree.dfLeaf2 <- tree.dfLeaf %>% filter(Drop==0) 
  
  tree.dfLeaf2 <- tree.dfLeaf2[,c("Lineage","name","Gene.OF")] 
  
  List.ReLin <- ggvita::ReLin(tree.dfLeaf2$Lineage,UseSubRoot = T)
  
  tree.dfLeaf2$Lineage <- List.ReLin$Leaf$NewLeafLin
  
  colnames(tree.dfLeaf2) <- c("Lineage","Name","Class")
  
  return(list("File.alm"=tree.dfLeaf2,"List.ReLin"=List.ReLin))
}
