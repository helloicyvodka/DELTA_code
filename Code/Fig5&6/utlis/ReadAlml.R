


ReadAlml<-function(file=""){
  
  
  
  the_result<-list()
  
  the_prefix<-c("Score","RootS","RootT","PruneS","PruneT","MatchS","MatchT","}","PValue","Min")
  the_text<-as.list(readLines(file))
  nl<-length(the_text)
  
  
  
  for(i in 1:nl){
    the_line<-the_text[[i]]
    if(regexpr("^[0-9]",the_line)){
      num2<-as.integer(regmatches(the_line,regexpr("^([0-9]+)",the_line)))
      if(length(num2)>0){
        num<-as.character(regmatches(the_line,regexpr("^([0-9]+)",the_line)))
        the_result[[num]]<-list()
        the_result[[num]]<-list("num"=c(num))
      }
    }
    
    
    for(i2 in 1:(length(the_prefix)-3)){
      if(startsWith(the_line,the_prefix[i2])){
        the_result[[num]][[as.character(the_prefix[i2])]]<-unlist(strsplit(the_line,split = ":"))[2]
      } 
    }
    if(startsWith(the_line,"PValue")){
      the_result[["PValue"]][["all"]]<-unlist(strsplit(the_line,split = ":"))[2]
    }
    if(startsWith(the_line,"Min")){
      the_result[["PValue"]][["Min"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[1],split = ":"))[2]
      the_result[["PValue"]][["Max"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[2],split = ":"))[2]
      the_result[["PValue"]][["AVG"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[3],split = ":"))[2]
    }
  }
  
  
  the_result[[1]]$Score
  
}