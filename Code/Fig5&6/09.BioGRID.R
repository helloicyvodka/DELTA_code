



BioGRID  <- 
  read.delim("./data/Biogrid/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.4.163.tab2.txt",sep="\t",header = T,stringsAsFactors = F)

BioGRID.2 <- BioGRID %>% filter(Experimental.System!="Two-hybrid"&Experimental.System.Type=="physical")

BioGRID.2 <- BioGRID.2 %>% mutate(Official.Symbol.Interactor.A = Official.Symbol.Interactor.A %>% toupper(),
                                  Official.Symbol.Interactor.B = Official.Symbol.Interactor.B %>% toupper())

BioGRID.2$Index.A <- sapply(1:nrow(BioGRID.2),function(x){
  
  r <- BioGRID.2[x,]
  My.all <- list(r$Systematic.Name.Interactor.A,r$Official.Symbol.Interactor.A,r$Synonyms.Interactor.A) %>% unlist() %>% as.character()
  My.all <- sapply(My.all,function(i)strsplit(i,split='\\|')) %>% unlist(recursive = T) %>% as.character()%>% gsub("CELE_","",x = .) 
  My.all[My.all!="-"] %>% unique() %>% toupper()
})

BioGRID.2$Index.B <- sapply(1:nrow(BioGRID.2),function(x){
  
  r <- BioGRID.2[x,]
  My.all <- list(r$Systematic.Name.Interactor.B,r$Official.Symbol.Interactor.B,r$Synonyms.Interactor.B) %>% unlist() %>% as.character()
  My.all <- sapply(My.all,function(i)strsplit(i,split='\\|')) %>% unlist(recursive = T) %>% as.character()%>% gsub("CELE_","",x = .) 
  My.all[My.all!="-"] %>% unique()%>% toupper()
})


Test.index <- sapply(1:nrow(BioGRID.2), function(x){
  
  r <- BioGRID.2[x,]
  
  m <- ifelse((intersect(unlist(r$Index.A),mut.genes) %>% length())!=0,1,0)
  n <- ifelse((intersect(unlist(r$Index.B),mut.genes) %>% length())!=0,1,0)
  
  m+n

})



BioGRID.3 <- BioGRID.2[Test.index==2,]

BioGRID.3$Index.A.2 <- sapply(BioGRID.3$Index.A,function(x){
  
  x[x %in% mut.genes]
  
})

BioGRID.3$Index.B.2 <- sapply(BioGRID.3$Index.B,function(x){
  
  x[x %in% mut.genes]
  
})

BioGRID.3 <- BioGRID.3 %>% filter(Index.A.2!=Index.B.2)



