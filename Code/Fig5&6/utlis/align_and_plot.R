library(devtools)
install_github("helloicyvodka/ggvita")  
library(ggvita);
library(plyr);
library(dplyr);
library(readr);
library(reshape2);
library(plyr);

##
## 1.prepare the data
##
## 1.1 load the dataset
## marker IDs as follow
## # RW10348 nhr-25 (hypodermal)
## # RW10425 pha-4 (pharynx)
## # RW10434 cnd-1 (neuronal)
## file format
## #Lineage cell_name average_expression binary_call time_point_expression
## #0 AB -802 0 -707,768,-485,-1431,-1479,-855,-1243,-2354,-827,-170,-906,65
## #00 ABa -1327 0 -4039,-605,-1229,-1446,-2319,-1935,-1531,-161,241,-464,-2058,-374
fMarker1.1 <- "../rawdata/digitalDev/WT/info_ZD_RW10348_WT_20110126_2_s2_emb1_edited.txt";
fMarker1.2 <- "../rawdata/digitalDev/WT/info_ZD_RW10425_WT_20100412_2_s1_emb1_edited.txt";
fMarker1.3 <- "../rawdata/digitalDev/WT/info_ZD_RW10434_WT_20110429_2_s1_emb1_edited.txt";
dfMarker1.1 <- read_delim(fMarker1.1,delim=" ",
                          col_names=T,col_types="ccdic");
dfMarker1.2 <- read_delim(fMarker1.2,delim=" ",
                          col_names=T,col_types="ccdic");
dfMarker1.3 <- read_delim(fMarker1.3,delim=" ",
                          col_names=T,col_types="ccdic");

fMarker2.1 <- "../rawdata/digitalDev/WT/info_ZD_RW10348_WT_20110126_2_s1_emb2_edited.txt";
fMarker2.2 <- "../rawdata/digitalDev/WT/info_ZD_RW10425_WT_20100412_2_s1_emb2_edited.txt";
fMarker2.3 <- "../rawdata/digitalDev/WT/info_ZD_RW10434_WT_20110429_2_s3_emb1_edited.txt";
dfMarker2.1 <- read_delim(fMarker2.1,delim=" ",
                          col_names=T,col_types="ccdic");
dfMarker2.2 <- read_delim(fMarker2.2,delim=" ",
                          col_names=T,col_types="ccdic");
dfMarker2.3 <- read_delim(fMarker2.3,delim=" ",
                          col_names=T,col_types="ccdic");

## 1.2 compile the data into lineage tree
## 
filterTerminal <- function(allLin) {
  useRows <- !( (paste(allLin,"0",sep="") %in% allLin) | allLin == "Root" );
  return(allLin[useRows]);
}
## 1.2 This section point to some issues that is corrected in 1.3, so just use 1.3
if(FALSE) { 
  
  filterTerminal <- function(dfLineage) {
    useRows <- !(paste(dfLineage$Lineage,"0",sep="") %in% dfLineage$Lineage);
    return(dfLineage[useRows,]);
  }
  ## 1.2.2 compile the tree
  dfMarker1 <- dfMarker1.1 %>% filterTerminal %>% dplyr::select(Lineage,nhr=binary_call) %>% 
    merge(dfMarker1.2 %>% filterTerminal %>% dplyr::select(Lineage,pha=binary_call)) %>% 
    merge(dfMarker1.3 %>% filterTerminal %>% dplyr::select(Lineage,cnd=binary_call));
  dfMarker2 <- dfMarker2.1 %>% filterTerminal %>% dplyr::select(Lineage,nhr=binary_call) %>% 
    merge(dfMarker2.2 %>% filterTerminal %>% dplyr::select(Lineage,pha=binary_call)) %>%
    merge(dfMarker2.3 %>% filterTerminal %>% dplyr::select(Lineage,cnd=binary_call));
  dfMarker1 %>% melt %>% 
    rbind.fill(dfMarker2 %>% melt %>% mutate(variable=paste(variable,"2",sep=""))) %>% 
    ggplot(aes(x=Lineage,y=variable,fill=factor(value))) + geom_tile()
}

## 1.3 Noticed missing nodes, they should inherit their mother cell's binary call
fillMissing1 <- function(dfLin,reqTerm){
  oriTerm <- setdiff(reqTerm,dfLin$Lineage);
  missTerm <- oriTerm;
  rowToFill <- rep(NA,length(missTerm));
  ## just resort to mother of terminal for now, don't trace further up
  #while(any(is.na(rowToFill))) {
    missTerm <- substr(missTerm,1,nchar(missTerm)-1),missTerm);
    rowToFill <- match(missTerm,dfLin$Lineage);
  #}
  data.frame(Lineage=oriTerm,binary_call=dfLin$binary_call[rowToFill]) %>%
    rbind.fill(dfLin) %>%
    return()
}

## 1.4 check whether the tree is complete
findIncompleteBranch <- function(allLin) {
  preNumLin <- 9999;
  while(preNumLin != length(allLin)) {
    preNumLin <- length(allLin);
    allMother <- allLin %>% substr(1,nchar(.)-1);
    freqTable <- allMother %>% table;
    allLin <- c(names(freqTable[freqTable == 2]),
                allLin[allMother %in% names(freqTable[freqTable == 1])])
  }
  return(allLin);
}
## 1.5 If not, make it complete by adding nodes
makeTreeComplete <- function(allLin) {
  incomBr <- allLin %>% findIncompleteBranch();
  while(incomBr %>% length() != 1) {
    theirDepth <- incomBr %>% nchar();
    maxDepth <- theirDepth %>% max();
    oriTail <- incomBr[theirDepth == maxDepth] %>% substr(nchar(.),nchar(.));
    oriHead <- incomBr[theirDepth == maxDepth] %>% substr(1,nchar(.)-1);
    addNode <- paste(oriHead,ifelse(oriTail == "0","1","0"),sep="");
    allLin <- c(allLin,addNode);
    incomBr <- allLin %>% findIncompleteBranch();
  }
  return(allLin)
}

## 1.6 finally, compile the tree
requiredTerminal <- unique(c(dfMarker1.1$Lineage,
                             dfMarker1.2$Lineage,
                             dfMarker1.3$Lineage,
                             dfMarker2.1$Lineage,
                             dfMarker2.2$Lineage,
                             dfMarker2.3$Lineage) ) %>% 
  filterTerminal %>%
  makeTreeComplete();
dfMarker1 <- dfMarker1.1 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,nhr=binary_call) %>% 
  merge(dfMarker1.2 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,pha=binary_call)) %>% 
  merge(dfMarker1.3 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,cnd=binary_call));
dfMarker2 <- dfMarker2.1 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,nhr=binary_call) %>% 
  merge(dfMarker2.2 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,pha=binary_call)) %>%
  merge(dfMarker2.3 %>% fillMissing1(requiredTerminal) %>% filter(Lineage %in% requiredTerminal) %>% dplyr::select(Lineage,cnd=binary_call));
dfMarker1 %>% melt %>% 
  rbind.fill(dfMarker2 %>% melt %>% mutate(variable=paste(variable,"2",sep=""))) %>% 
  ggplot(aes(x=Lineage,y=variable,fill=factor(value))) + geom_tile();

## ?? 1.7 
## The nodes that is added to make the tree complete is assumed its mother's expression state
## Need more statistical correction for the expression ??

##
## 2.do the alignment
##
## # 2.1 prepare the score matrix
ndiff <- function(x,y){
  sx <- strsplit(x,"") %>% unlist();
  sy <- strsplit(y,"") %>% unlist();
  if(length(sx) != length(sy)) {return(NA);}
  return(length(which(sx != sy)));
}
dfTypes <- dfMarker1 %>% 
  rbind.fill(dfMarker2) %>%
  group_by(nhr,pha,cnd) %>%
  sample_n(1) %>%
  mutate(type=paste(nhr,pha,cnd,sep="")) %>%
  ungroup() %>%
  dplyr::select(-Lineage);
dfScore <- dfTypes$type %>% 
  combn(2) %>%
  t() %>%
  data.frame(stringsAsFactors = F) %>%
  rbind.fill(data.frame(X1=dfTypes$type,X2=dfTypes$type)) %>%
  group_by(X1,X2) %>%
  mutate(score = -exp(ndiff(X1,X2))) %>%
  ungroup() %>%
  mutate(score = round((score - mean(score,na.rm=T)) / sd(score,na.rm=T) * 10));
dfScore$score[grep("NA",dfScore$X1)] <- 0;
dfScore$score[grep("NA",dfScore$X2)] <- 0;

dfScore[dfScore$X1 != dfScore$X2,"score"] <- -100;
dfScore[dfScore$X1 == dfScore$X2,"score"] <- 10;

fnCost <- "01.score.tsv";
write_tsv(dfScore,fnCost,col_names=F);

## 2.2 prepare the tree
dfTree1 <- dfMarker1 %>% 
  ungroup %>% 
  mutate(type=paste("E",nhr,pha,cnd,sep="")) %>% 
  mutate(Name = sample.int(nrow(.)) ) %>%
  dplyr::select(Lineage,Name,Class=type)
fnTree1 <- "01.treeS.alm";
write_tsv(dfTree1,
          fnTree1);
dfTree2 <- dfMarker2 %>% 
  ungroup %>% 
  mutate(type=paste("E",nhr,pha,cnd,sep="")) %>% 
  mutate(Name = sample.int(nrow(.)) ) %>%
  dplyr::select(Lineage,Name,Class=type)
fnTree2 <- "01.treeT.alm";
write_tsv(dfTree2,
          fnTree2);
fnRes <- "01.res.txt";
system(sprintf("perl -pi -e 'chomp if eof' %s",fnTree1)); ## remove the last "newline"
system(sprintf("perl -pi -e 'chomp if eof' %s",fnTree2)); ## otherwise DELTA will fail
system(sprintf("perl -pi -e 'chomp if eof' %s",fnCost));
system(sprintf("../code/DELTA/bin/Release/DELTA -treeS %s -treeT %s -cost %s -method l -prune 100 -max_target 100 -test 0 -outfile %s",
               fnTree1,fnTree2,fnCost,fnRes));

##
## 3.plot one of the symmetric subtree pairs
##
readal.alm(fnTree1,fnTree2);
res.Delta <- readal.alml2(fnRes);
p <- ggvita(res.Delta$`2`,print=T);
allTypes <- c(dfTree1$Class,dfTree2$Class) %>% unique;
type2col <- brewer.pal(length(allTypes),"Set2");
names(type2col) <- allTypes;
p + 
  geom_tippoint(aes(color=label)) + 
  scale_color_manual(values=type2col)

##
## 4. check the rank of the symmtry pairs
##
## 4.1 load the sym pair list
dfSym <- read_csv("../rawdata/Cel.symmetry/symmetry-sisters.csv",col_types="cccc");
dfDelta.root <- lapply(res.Delta,function(x){data.frame(s.root=x$RootS,t.root=x$RootT)}) %>% rbind.fill;

match1.s <- match(dfDelta.root$s.root,dfSym$Sister1_binary);
match1.t <- match(dfDelta.root$t.root,dfSym$Sister2_binary);
match1 <- match1.s[match1.s == match1.t];
match2.s <- match(dfDelta.root$s.root,dfSym$Sister2_binary);
match2.t <- match(dfDelta.root$t.root,dfSym$Sister1_binary);
match2 <- match2.s[match2.s == match2.t];
foundSister <- ifelse(is.na(match1),match2,match1);
lapply(res.Delta,function(x){x$Score}) %>% unlist %>% as.numeric();

lapply(paste("^",dfSym$Sister1_binary,sep=""),grep,dfTree1$Lineage,perl=T);

if(FALSE) {
  readal.alml2<-function(file){
    
    if(!exists("alm_label")){
      stop("Please run readal.alm firstly and create the alm_label varaiable! ATTENTION: readal.alm will creates alm_label variable automately.")
    }
    
    
    the_result<-list()
    
    the_prefix<-c("Score","RootS","RootT","PruneS","PruneT","MatchS","MatchT","}","PValue","Min")
    
    the_text<-readLines(file)
    
    nl<-length(the_text)
    
    
    
    for(i in 1:(nl-1)){
      
      the_line<-the_text[[i]]
      
      
      if(regexpr("^[0-9]",the_line)==T){
        
        num<-as.character(regmatches(the_line,regexpr("^([0-9]+)",the_line)))
        
        if(startsWith(the_text[[(i+1)]],"Score")==T){
          
          if(length(num)>0){
            the_result[[num]]<-list()
            the_result[[num]]<-list("score_order"=as.numeric(num))
          }
        }
        
      }else{
        
        for(i2 in 1:(length(the_prefix)-3)){
          
          if(startsWith(the_line,the_prefix[i2])){
            
            the_result[[num]][[as.character(the_prefix[i2])]]<-unlist(strsplit(the_line,split = ":"))[2]
            
          }
        }
        
        
        if(startsWith(the_line,"PValue")){
          
          the_result[[num]][["PValue"]]<-list()
          
          the_result[[num]][["PValue"]][["all"]]<-unlist(strsplit(the_line,split = ":"))[2]
          
        }
        
        if(startsWith(the_line,"Min")){
          
          the_result[[num]][["PValue"]][["Min"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[1],split = ":"))[2]
          
          the_result[[num]][["PValue"]][["Max"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[2],split = ":"))[2]
          
          the_result[[num]][["PValue"]][["AVG"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[3],split = ":"))[2]
          
          the_result[[num]][["PValue"]][["pvalue"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[4],split = ":"))[2]
          
        }
      }
      
    }
    
    
    for(i in 1:length(the_result)){
      class(the_result[[i]])<-c("alml",class(the_result[[i]]))
    }
    
    class(the_result)<-c("alml_list",class(the_result))
    
    return(the_result)
    
  }
}
