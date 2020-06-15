
# Content
# . Find.root
# . Find.tips
# . Find.missed.mother
# . Find.missed.sister
# . Find.addNode
# . Make.tree.complete
# . Find.unpaired.tip
# . LN_to_Bin.2 ( "Root" --> "" )


library(plyr)
library(dplyr)
library(testit)

#----------------------------------------------------------------------------

Find.root <- function(allLin){
  
  the.ava.roots <- c("Root","root","")
  
  if(allLin[allLin %in% the.ava.roots] %>% length()  >1){
    
    stop("The lineage sequence includes abnormal roots!")
  }
  
  if("Root" %in% allLin){
    
    the.root <- "Root"
    
  }else if( "root" %in% allLin){
    
    the.root <- "root"
    
  }else{
    
    the.min.depth <- allLin %>% nchar() %>% min()
    the.potential.root <- allLin[nchar(allLin)==the.min.depth]
    the.potential.root.son <- allLin[nchar(allLin)==the.min.depth+1]
    
    if(length(the.potential.root)==2){
      the.root <- the.potential.root %>% unique() %>% substr(1,nchar(.)-1) %>% unique()
      warning("This lineage sequence doesn't include the root lineage!")
    }else if(length(the.potential.root)==1){
      if(length(the.potential.root.son)==2){
        the.root <- the.potential.root
      }else{stop("The tree doesn't have a normal root!")}
    }else{stop("The tree doesn't have a normal root!")}
    
    
  }
  
  
  return(the.root)
  
}



#----------------------------------------------------------------------------
Find.tips <- function(allLin){
  
    
    the.root <- Find.root(allLin)
    
    if(has_warning(Find.root(allLin))==T){
      
      allLin.2 <- allLin 
      
    }else{
      
      allLin.2 <- allLin[allLin!=the.root]
      
    }
    

  
 
  the.allMother <- allLin.2 %>% substr(1,nchar(.)-1) %>% unique()
  the.tips <- setdiff(allLin.2,the.allMother)
  
  return(the.tips)
  
}


#----------------------------------------------------------------------------

Find.missed.mother <- function(allLin){
  
  the.root <- Find.root(allLin)
  
  if(has_warning(Find.root(allLin))==T){
    
    allLin.2 <- allLin 
    
  }else{
    
    allLin.2 <- allLin[allLin!=the.root]
    
  }
  
  all.missed <- c()
  
  
  repeat{
    
    
    
    allMother.2 <- allLin.2 %>% substr(1,nchar(.)-1)
    
    the.missed.mother <- setdiff(allMother.2,allLin.2)
    the.missed.mother <- setdiff(the.missed.mother,
                                 ifelse(the.root %in% c("Root","root"),"",the.root))
    
    if(length(the.missed.mother)==0)break
    
    allLin.2 <- c(allLin.2,the.missed.mother) %>% unique()
    all.missed <- c(all.missed,the.missed.mother) %>% unique()
    
  }
  
  
  return(all.missed)
}



#----------------------------------------------------------------------------


Find.missed.sister <- function(allLin){
  
  the.root <- Find.root(allLin)
  
  if(has_warning(Find.root(allLin))==T){
    
    allLin.2 <- allLin 
    
  }else{
    
    allLin.2 <- allLin[allLin!=the.root]
    
  }
  
   
  oriTail <- allLin.2 %>% substr(nchar(.),nchar(.));
  oriHead <- allLin.2 %>% substr(1,nchar(.)-1);
  allLin.2.mirror <- paste(oriHead,ifelse(oriTail == "0","1","0"),sep="");
  
  the.missed.sister <- setdiff(allLin.2.mirror,allLin.2)
  
  return(the.missed.sister)
}


#----------------------------------------------------------------------------



Find.addNode <- function(allLin){
  
  addNode <- c()
  
  repeat{
    
    incomBr.1 <- allLin %>% Find.missed.sister();
    incomBr.2 <- allLin %>% Find.missed.mother();
    incomBr <- unique(c(incomBr.1,incomBr.2));
    
    if(incomBr %>% length()== 0) break
    
    addNode <- c(addNode,incomBr)
    allLin <- unique(c(allLin,incomBr));
    
    
  }
  
  return(addNode)
  
}

#----------------------------------------------------------------------------


Find.unpaired.tip <- function(allLin){
    
  the.root <- Find.root(allLin)
  
  if(has_warning(Find.root(allLin))==T){
    
    allLin.2 <- allLin 
    
  }else{
    
    allLin.2 <- allLin[allLin!=the.root]
    
  }
  
  addNode <- 
    c(allLin,allLin %>% Find.missed.mother()) %>% 
    unique() %>% 
    Find.missed.sister()
  
  
  if(length(addNode)!=0){
    

    oriTail <- addNode %>% substr(nchar(.),nchar(.));
    oriHead <- addNode %>% substr(1,nchar(.)-1);
    the.unpaired.tip <- paste(oriHead,ifelse(oriTail == "0","1","0"),sep="");

    return(the.unpaired.tip)
  }else{
    
    warning("There is no unpaired tip")
    
  }
  
  
  
}

#----------------------------------------------------------------------------


Make.tree.complete <- function(allLin) {
  
  allLin <- c(allLin,Find.addNode(allLin))
  return(allLin)
  
}


#----------------------------------------------------------------------------
Find.tree.from.tips <- function(allTips){
  
  the.above.nodes <- allTips
  repeat{
    
    
    if(the.above.nodes %>% length()>1){
      
      the.above.nodes <- the.above.nodes %>% substr(1,nchar(.)-1) %>% unique()
      
      allTips <- c(allTips,the.above.nodes) %>% unique()
      
    }else{break}
    
    
    
  }
  allTips
}


#----------------------------------------------------------------------------



Na.omit.expr.dt <- function(com.expr){
  

  repeat{
    
    com.expr.tips <- com.expr %>% filter(Lineage %in% (Lineage %>% Find.tips())) %>% na.omit() %>% `$`("Lineage")
    
    com.expr.tips.na <- setdiff(com.expr %>% filter(Lineage %in% (Lineage %>% Find.tips())) %>% `$`("Lineage"),
                                com.expr.tips)
    
    if(com.expr.tips.na %>% length()!=0){
      
      com.expr <- com.expr  %>% filter(!Lineage %in% com.expr.tips.na)
      
    }else{break}
    
  }
  
  
  if(com.expr$Lineage %>% Find.addNode() %>% length()!=0)stop("Error: Tree is not complete!")
  
  return(com.expr)
  
}
#----------------------------------------------------------------------------


LN_to_Bin.2 <- function(the_LN_vec) {
  
  
  
  #########
  
  the_prefix <- c(
    "AB",
    "Ab",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G1",
    "G2",
    "H1L",
    "H1R",
    "H2L",
    "H2R",
    "K",
    "M",
    "MS",
    "P0",
    "P1",
    "P10",
    "P11",
    "P12",
    "P2",
    "P3",
    "P4",
    "P5",
    "P6",
    "P7",
    "P8",
    "P9",
    "QL",
    "QR",
    "TL",
    "TR",
    "U",
    "V1L",
    "V1R",
    "V2L",
    "V2R",
    "V3L",
    "V3R",
    "V4L",
    "V4R",
    "V5L",
    "V5R",
    "V6L",
    "V6R",
    "W",
    "Y",
    "Z",
    "Z1",
    "Z4",
    "Z2",
    "Z3",
    "EMS"
  )
  
  #########
  
  the_true_prefix <- c(
    "Za", #"AB"
    "Za",
    "Zaprppppapa",
    "Zppa",
    "Zpppa",
    "Zpap",
    "Zaplppppapp",
    "Zaprpaaaapa",
    "Zaplapaapa",
    "Zaplaaappp",
    "Zaarpapppp",
    "Zaarppaaap",
    "Zaarpppaap",
    "Zaplpapppaa",
    "Zpaaapaapp",
    "Zpaa",
    "Z",
    "Zp",#P1(CXL:"Zaplapaapp"->YM:"Zp")
    "Zaprapapap", #P10
    "Zaplapappa", #P11
    "Zaprapappa", #P12
    "Zpp",#P2(CXL:"Zaprapaapp"->YM:"Zpp")
    "Zppp", #P3(CXL: "Zaplappaaa"->YM:"Zppp")
    "Zpppp",#P4(CXL:"Zaprappaaa"->YM:"Zpppp")
    "Zaplappaap",
    "Zaprappaap",
    "Zaplappapp",
    "Zaprappapp",
    "Zaplapapap",
    "Zaplapapaaa",
    "Zaprapapaaa",
    "Zaplappppp",
    "Zaprappppp",
    "Zaplppppapa",
    "Zaarppapaa",
    "Zaarppppaa",
    "Zaarppapap",
    "Zaarppppap",
    "Zaplappapa",
    "Zaprappapa",
    "Zaarppappa",
    "Zaarpppppa",
    "Zaplapapaap",
    "Zaprapapaap",
    "Zaarppappp",
    "Zaarpppppp",
    "Zaprapaapa",
    "Zaprpppaaaa",
    "Z",
    "Zpaapppaap",
    "Zpaaappaap",
    "Zppppp",
    "Zppppa",
    "Zpa"
  )
  
  
  
  the_prefix_df <-
    data.frame(Prefix = the_prefix, True_prefix = the_true_prefix) %>%
    list.parse()
  
  
  ###########
  
  
  the_LN_2_trueLN <- function(x) {
    
    the_matched_prefix <- the_prefix[startsWith(x, the_prefix)]
    
    the_matched_prefix <-
      the_matched_prefix[which(nchar(the_matched_prefix) == (the_matched_prefix %>% nchar() %>% max()))]
    
    if (length(the_matched_prefix) == 0) {
      stop("The prefix was not matched!")
    }
    
    #a -> 0; p -> 1
    #l -> 0; r -> 1
    #d -> 0; v -> 1
    
    x <- x %>%
      as.character() %>%
      gsub("\\.", "", .) %>%
      gsub(" ", "", .) %>%
      sub(
        the_matched_prefix,
        the_prefix_df %>>%
          list.filter(Prefix == the_matched_prefix) %>>%
          list.mapv(True_prefix)
        ,
        .
      )
    
    if(x!="Z"){
      
      x<-sub("Z", "", x) %>%
        gsub("a", "0", .) %>%
        gsub("l", "0", .) %>%
        gsub("d", "0", .) %>%
        gsub("p", "1", .) %>%
        gsub("r", "1", .) %>%
        gsub("v", "1", .)
      
    }else{
      x<-""
    }
    
    return(x)
  }
  
  ###########
  
  
  return(mclapply(the_LN_vec,
                  the_LN_2_trueLN,
                  mc.cores = 16) %>%
           unlist())
  
  
}


















