library(plyr)
library(dplyr)


find.Incomplete.Branch <- function(allLin) {
  
  repeat{
    
    preNumLin <- length(allLin);
    
    allMother <- allLin %>% substr(1,nchar(.)-1);
    
    freqTable <- allMother %>% table;
    
    allLin <- c(names(freqTable[freqTable == 2]),
                
                allLin[allMother %in% names(freqTable[freqTable == 1])])
    
    if(preNumLin == length(allLin))break
  }
  
  return(allLin[! allLin %in% c("Root","")]);
  
}

find.Incomplete.Branch.level <- function(allLin,level) {
  
  for(i in 1:level+1){
    
    preNumLin <- length(allLin);
    
    allMother <- allLin %>% substr(1,nchar(.)-1);
    
    freqTable <- allMother %>% table;
    
    allLin <- c(names(freqTable[freqTable == 2]),
                
                allLin[allMother %in% names(freqTable[freqTable == 1])])
    
    
  }
    
    
  
  
  return(allLin[! allLin %in% c("Root","")]);
  
}
