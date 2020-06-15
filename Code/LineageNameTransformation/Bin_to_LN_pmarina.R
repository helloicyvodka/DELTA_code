
library(xlsx)
library(rlist)
library(pipeR)
library(testit)
library(magrittr)

PP <- data.table::fread("~/2017-2018/tree/Pmarina_Celegans/data/pma-new.alm",colClasses = "character") 
colnames(PP) <- c("Bin","LN_clean","CellType")

pp_list <- list.parse(PP[,c("LN_clean","Bin")])

Bin_to_LN_pma <- function(i){
  
  if(is.na(i))return(NA)
  if(i=="")return(NA)
  
  i <- as.character(i)
  if(i=="1")return("EMS")

  t <- 1  
  suffix <- ""
  
  repeat{
    
    ilist<- pp_list %>>% 
      list.filter(Bin %>% as.character %>% startsWith(prefix=i))
    if(length(ilist)!=0){break}
    i_s <- substr(i,nchar(i)-1,nchar(i))
    suffix <- paste0(suffix,ifelse(i_s=="0","l","r"))
    i <- substr(i,1,nchar(i)-1)
    
  }
  
  
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

