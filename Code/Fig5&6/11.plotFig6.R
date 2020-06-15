library(dplyr)
library(reshape2)
library(ggplot2)

midDELTA <- 5005;

nBoot <- 1000;

x=0

toPlot.ratio <- lapply(c(0:8)*300,function(x){
  
  topMean <- toPlot.worm.2 %>%
    
    dplyr::filter(DELTA.score > (midDELTA + x) ) %>%
    
    dplyr::select(dNdS,dN,dS,logexpr) %>%
    
    colMeans(na.rm=T);
  
  
  
  bottomMean <- toPlot.worm.2 %>%
    
    filter(DELTA.score < (midDELTA - x) ) %>%
    
    select(dNdS,dN,dS,logexpr) %>%
    
    colMeans(na.rm=T);
  
  
  
  myRes <- data.frame(idx = sample.int(nrow(toPlot.worm.2),nrow(toPlot.worm.2) * nBoot,replace=T),
                      
                      bootGrp = c(1:nBoot)) %>%
    
    dplyr::group_by(bootGrp) %>%
    
    do({
      
      myDf <- .;
      
      topMean2 <- toPlot.worm.2[myDf$idx,] %>%
        
        filter(DELTA.score > (midDELTA + x) ) %>%
        
        select(dNdS,dN,dS,logexpr) %>%
        
        colMeans(na.rm=T);
      
      bottomMean2 <- toPlot.worm.2[myDf$idx,] %>%
        
        filter(DELTA.score < (midDELTA - x) ) %>%
        
        select(dNdS,dN,dS,logexpr) %>%
        
        colMeans(na.rm=T);
      
      data.frame(top=topMean2,bottom=bottomMean2) %>%
        
        t() %>%
        
        melt() %>%
        
        rename(`group`=`Var1`,`var`=`Var2`)
      
    }) %>%
    
    group_by(group,var) %>%
    
    dplyr::summarise(sd = sd(value,na.rm=T),
                     
                     se = sd(value,na.rm=T)/sqrt(1000),
                     
                     ciup = quantile(value,prob=c(0.975),na.rm=T),
                     
                     cidown = quantile(value,prob=c(0.0255),na.rm=T));
  
  myRes %>%
    
    ungroup() %>%
    
    mutate(obs = unlist(c(topMean,bottomMean)),
           
           x = x,
           
           topCut = midDELTA + x,
           
           bottomCut = midDELTA - x) %>%
    
    mutate(label = ifelse(group == "top", sprintf(">%.0f",topCut),sprintf("<%.0f",bottomCut))) %>%
    
    return()
  
}) %>% plyr::rbind.fill();



toPlot.ratio %>%
  
  ggplot(aes(x=x,y=obs,group=group,fill=group)) +
  
  geom_bar(stat="identity",color="black",position=position_dodge(),width=150) +
  
  geom_errorbar(aes(ymax=ciup,ymin=cidown),width=100,position=position_dodge(width=150)) +
  
  facet_grid(var ~ .,scales = "free",) +
  
  scale_x_continuous("Groups of DELTA score") +
  
  scale_y_continuous("")+
  
  theme_classic() +
  
  #style.print() +
  
  theme(axis.text.x = element_text(angle = 80,hjust=.5,vjust=0.5))
