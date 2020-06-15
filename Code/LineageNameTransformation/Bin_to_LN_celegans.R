
library(xlsx)
library(rlist)
library(pipeR)
library(testit)
library(magrittr)

CE <- data.table::fread("~/2017-2018/tree/Pmarina_Celegans/data/fun.alm",colClasses = "character") 
colnames(CE) <- c("Bin","LN_clean","CellType")

ce_list <- list.parse(CE[,c("LN_clean","Bin")])

Bin_to_LN_cel <- function(i){
  
  if(is.na(i))return(NA)
  if(i=="")return(NA)
  
  i <- as.character(i)

  t <- 1  
  
  ilist<- ce_list %>>% 
      list.filter(Bin %>% as.character %>% startsWith(prefix=i)) 
  
  repeat{
    

    
    ilist_s <- ilist[t]
    
    ii <- ilist_s %>>%
      list.select(Bin) %>>%
      unlist() %>>% 
      as.character()
    
    iLN <- ilist_s %>>%
      list.select(LN_clean)%>>%
      unlist() %>>% 
      as.character()
    
    ttt <-nchar(iLN)-(nchar(ii)-nchar(i))
    t <- t+1
    if(ttt>0){break}
    
  }
  
  return(substr(iLN,1,ttt))
}

