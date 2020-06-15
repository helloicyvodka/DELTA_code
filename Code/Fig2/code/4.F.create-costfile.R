library(dplyr)

F.create.costfile <- function(tree1.df,score,m=1,u=1,silent=FALSE){
  
  class.uni <- tree1.df$Class %>% unique()
  
  dfscore <- expand.grid(class.uni,class.uni) %>% data.frame(stringsAsFactors = F)
  
  colnames(dfscore) <- c("S","T")
  
  if (score == "common"){
    # Hamming distance 
    dfscore$Score <- apply(dfscore,1,function(x){
      hd <- sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],"")));
      m*hd-u*(8-hd)
    })
    
  }else if (score == "same"){
    dfscore$Score <- (dfscore$S == dfscore$"T")
    dfscore$Score <-  ifelse(dfscore,m,-u)

  }else if (score == "norm"){
    dfscore$Score <- apply(dfscore,1,function(x) sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],""))))
    dfscore$Score <- round((dfscore$Score - mean(dfscore$Score))/sd(dfscore$Score)) + 1
  }else if (!silent){print("invalid value for score, score can be same, common or norm")}
  
  
  as.matrix(dfscore)
  
  
}
  