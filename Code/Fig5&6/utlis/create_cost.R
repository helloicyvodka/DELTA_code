

create.cost <- function(mt.s, unmt.s,file){
  
  
  tr.Class <- c(0:7)
  
  cost.mt <-as.data.frame(expand.grid(tr.Class,tr.Class))
  
  
  
  
  
  Dec2Bin <- function(x) {
    i <- 0
    string <- numeric(32)
    if(x==0){return(x)}
    while(x > 0) {
      string[32 - i] <- x %% 2
      x <- x %/% 2
      i <- i + 1
    }
    first <- match(1, string)
    paste(as.character(string[first:32]),collapse="")
  }
  
  
  
  to.bin.vec <- function(x){
    
    x <- Dec2Bin(x)
    
    if(nchar(x) == 1){
      x <- as.vector(as.numeric(x))
    }else{
      x <- as.numeric(unlist(strsplit(x, split = "")))
    }
    
    
    if( length(x) <3 ){
      x <- c(rep(0,3-length(x)),x)
    }
    
    return(x)
  }
  
  
  
  
  
  
  cost.mt$score <- sapply(1:nrow(cost.mt),function(x){
    
    r <- cost.mt[x,]
    
    r.1 <- to.bin.vec(r[1])
    r.2 <- to.bin.vec(r[2])
    
    r.r <- abs(r.1-r.2)
    
    y <- length(r.r[r.r==0])*mt.s+length(r.r[r.r==1])*(-unmt.s)
    
    return(y)
    
  })
  
  write.table(as.matrix(cost.mt),file = file ,quote = F,row.names = F,col.names = F)
  
  
}


write.table(as.matrix(cost.mt),file = "./data/data_from_digital/cost_mt10_unmt1.tsv",quote = F,row.names = F,col.names = F)




