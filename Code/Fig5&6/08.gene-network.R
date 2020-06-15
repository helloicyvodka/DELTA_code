
library(dplyr)
library(ggvita)
library(ggplot2)
library(cowplot)

String.ineraction <- read.table("./data/string/string_interactions.tsv",header = T,stringsAsFactors = F)
String.ineraction <- String.ineraction %>% mutate(node1=toupper(node1),node2=toupper(node2))



mut.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/mut-tree-tips",full.names = T)

#wt.files.list <- list.files("~/2017-2018/tree/tree_fate_change/data2/wt-tree-tips",full.names = T)

cost <- "/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data/data_from_digital/cost2018.tsv"

#String.outfile.parent.folder <- "./data3/String"

# result.list.files.genes  <- unlist(sapply(list.files("./data3/same-almg/mut-vs-mut/",full.names = F),
#                                         function(x){strsplit(x,split = "_")[[1]][1]})) %>% as.character()







String.ineraction <- String.ineraction %>% filter(node1 %in% mut.genes & node2 %in% mut.genes)

library(parallel)


All.mut <- 
  expand.grid(mut.genes,
              mut.genes,eins
              stringsAsFactors = F) %>% 
  dplyr::rename(`node1`=`Var1`,`node2`=`Var2`) %>% 
  dplyr::filter(node1>node2)





mclapply(1:nrow(All.mut),function(i){
  
  r <- All.mut[i,]
  
  output <- paste0("./data3/gene_network_pr1_all/",i)
  
  dir.create(output)
  
  mut.A <- mut.full.df.list[[r$node1]]
  
  mut.B <- mut.full.df.list[[r$node2]]
  
  if(is.null(mut.B)|is.null(mut.A)){
    
  }else{
    
    AB.inter <- intersect(mut.A$Lineage,mut.B$Lineage)
    
    AB.inter <- c(c("Root","0","1","00","01","11","10"),AB.inter) %>% unique()
    
    if(is.null(ggvita::Find.addNode(AB.inter))){
      
      AB.tips <- ggvita::Find.tips(AB.inter)
      
      mut.A.2 <-data.frame(Lineage=mut.A$Lineage,Name=mut.A$cell_name,Class=mut.A$Class.num,stringsAsFactors = F)
      mut.B.2 <-data.frame(Lineage=mut.B$Lineage,Name=mut.B$cell_name,Class=mut.B$Class.num,stringsAsFactors = F)
      
      ggvita::write.alm(mut.A.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node1,".alm"))
      ggvita::write.alm(mut.B.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node2,".alm"))
      
      ggvita::DELTA(
        treeS = paste0(output,"/",r$node1,".alm"),
        treeT = paste0(output,"/",r$node2,".alm"),
        cost = cost,
        outfile = paste0(output,"/","My.almg"),
        method = "g",
        prune=1,
        DELTA.address =NULL,
        test = 100
      )
    }
  }
},mc.cores = 30)
    
      



All.mut$DELTA.score <- 
  mclapply(1:nrow(All.mut),function(i){
    
    r <- All.mut[i,]
    
    output <- paste0("./data3/gene_network_pr1_all/",i)
    
    #dir.create(output)
    
    mut.A <- mut.full.df.list[[r$node1]]
    mut.B <- mut.full.df.list[[r$node2]]
    
    if(is.null(mut.B)|is.null(mut.A)){
      
      return(NA)
      
    }else{
      
    
        t <- 
          ggvita::readal(
            fileS = paste0(output,"/",r$node1,".alm"),
            fileT = paste0(output,"/",r$node2,".alm"),
            outfile = paste0("./data3/gene_network_pr1_all/",i,"/My.almg"),
            cost=cost,
            method = "g"
          )
        

        return(ifelse(is.null(t$Score),NA,t$Score) )     
               
        }
    }
    
    
  ,mc.cores=30) 



All.mut$DELTA.score <-  
  All.mut$DELTA.score %>% unlist()








String.result.common <-

parallel::mclapply(1:nrow(String.ineraction),FUN = function(i){

    r <- String.ineraction[i,]


    output <- paste0("./data3/gene_network_pr40/",i)

    #dir.create(output)


    mut.A <- mut.full.df.list[[r$node1]]
    mut.B <- mut.full.df.list[[r$node2]]

    if(is.null(mut.B)|is.null(mut.A)){

      return(NA)

    }else{

      AB.inter <- intersect(mut.A$Lineage,mut.B$Lineage)

      AB.inter <- c(c("Root","0","1","00","01","11","10"),AB.inter) %>% unique()

      if(is.null(ggvita::Find.addNode(AB.inter))){

        AB.tips <- ggvita::Find.tips(AB.inter)

        mut.A.2 <-data.frame(Lineage=mut.A$Lineage,Name=mut.A$cell_name,Class=mut.A$Class.num,stringsAsFactors = F)
        mut.B.2 <-data.frame(Lineage=mut.B$Lineage,Name=mut.B$cell_name,Class=mut.B$Class.num,stringsAsFactors = F)

        ggvita::write.alm(mut.A.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node1,".alm"))
        ggvita::write.alm(mut.B.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node2,".alm"))

        ggvita::DELTA(
          treeS = paste0(output,"/",r$node1,".alm"),
          treeT = paste0(output,"/",r$node2,".alm"),
          cost = cost,
          outfile = paste0(output,"/","My.almg"),
          method = "g",
          prune=1000,
          DELTA.address =NULL,
          test = 100
        )
        
        

        t <- 
        ggvita::readal(
          fileS = paste0(output,"/",r$node1,".alm"),
          fileT = paste0(output,"/",r$node2,".alm"),
          outfile = paste0("./data3/gene_network_pr40/",i,"/My.almg"),
          cost=cost,
          method = "g"
        )

        #t <- ggvita::readal.almg(paste0(output,"/","My.almg"))

        return(c(t)
      }
      }
    },mc.cores = 30,mc.silent = T)



for(i in 1:nrow(String.ineraction)){
  
  r <- String.ineraction[i,]
  
  output <- paste0("./data3/gene_network_pr1000/",i)
  
  dir.create(output)
  
  mut.A <- mut.full.df.list[[r$node1]]
  mut.B <- mut.full.df.list[[r$node2]]
  
  if(is.null(mut.B)|is.null(mut.A)){
    
    next()
    
  }else{
    
    AB.inter <- intersect(mut.A$Lineage,mut.B$Lineage)
    
    AB.inter <- c(c("Root","0","1","00","01","11","10"),AB.inter) %>% unique()
    
    if(is.null(ggvita::Find.addNode(AB.inter))){
      
      AB.tips <- ggvita::Find.tips(AB.inter)
      
      mut.A.2 <-data.frame(Lineage=mut.A$Lineage,Name=mut.A$cell_name,Class=mut.A$Class.num,stringsAsFactors = F)
      mut.B.2 <-data.frame(Lineage=mut.B$Lineage,Name=mut.B$cell_name,Class=mut.B$Class.num,stringsAsFactors = F)
      
      ggvita::write.alm(mut.A.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node1,".alm"))
      ggvita::write.alm(mut.B.2 %>% filter(Lineage %in% AB.tips),paste0(output,"/",r$node2,".alm"))
      
      ggvita::DELTA(
        treeS = paste0(output,"/",r$node1,".alm"),
        treeT = paste0(output,"/",r$node2,".alm"),
        cost = cost,
        outfile = paste0(output,"/","My.almg"),
        method = "g",
        prune=1000,
        DELTA.address =NULL,
        test = 100
      )
      
      t <- 
        ggvita::readal(
          fileS = paste0(output,"/",r$node1,".alm"),
          fileT = paste0(output,"/",r$node2,".alm"),
          outfile = paste0("./data3/gene_network_pr1/",i,"/My.almg"),
          cost=cost,
          method = "g"
        )
      
      String.result.common[[i]] <- t
    }
  }
  
   cat(i,"\n")  
}






String.ineraction$DELTA.score <- sapply(1:nrow(String.ineraction),function(i){
  
  # r <- String.ineraction[i,]
  # 
  # treeS <- mut.files.list[grep(as.character(r$node1), mut.files.list)]
  # treeT <- mut.files.list[grep(as.character(r$node2), mut.files.list)]
  # outfile.mut_mut <- paste0(as.character(r$node1),
  #                           "_",
  #                           as.character(r$node2),
  #                           "_.almg")
  # if(outfile.mut_mut %>% file.exists())return(NA)
  
  s <- String.result.common[[i]]$Score
  
  ifelse(is.null(s),NA,s)
  
}) 


String.ineraction$DELTA.pvalue <- sapply(1:nrow(String.ineraction),function(i){
  
  # r <- String.ineraction[i,]
  # 
  # treeS <- mut.files.list[grep(as.character(r$node1), mut.files.list)]
  # treeT <- mut.files.list[grep(as.character(r$node2), mut.files.list)]
  # outfile.mut_mut <- paste0(as.character(r$node1),
  #                           "_",
  #                           as.character(r$node2),
  #                           "_.almg")
  # if(outfile.mut_mut %>% file.exists())return(NA)
  
  s <- String.result.common[[i]]$PValue$pvalue,$Score
  
  ifelse(is.null(s),NA,s)
  
}) 





View(String.ineraction) 




# String.ineraction.2 <- 
#   String.ineraction %>%
#   filter(!is.na(DELTA.score)) %>% 
#   mutate(neg.log10.pvalue=-log10(as.numeric(DELTA.pvalue)))


String.ineraction.3 <- 
  String.ineraction %>%
  filter(!is.na(DELTA.pvalue)&DELTA.pvalue>0) %>% 
  mutate(neg.log10.pvalue=-log10(as.numeric(DELTA.pvalue)))



# "node1"                                 "node2"                                 "node1_string_internal_id"              "node2_string_internal_id"             
# "node1_external_id"                     "node2_external_id"                     "neighborhood_on_chromosome"            "gene_fusion" 

                                                                   
# "phylogenetic_cooccurrence"                        
#"homology"                              
#"coexpression"                          
#"experimentally_determined_interaction"
# "database_annotated"                    
#"automated_textmining"                  
#"combined_score"                        
#"DELTA.score"                          
#"DELTA.pvalue"                          
#"Sub.group"                          


# phylogenetic_cooccurrence 

cor.test(String.ineraction.2$DELTA.score %>% as.numeric(),
         String.ineraction.2$phylogenetic_cooccurrence %>% as.numeric(),
         method = "sp")

"No"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$phylogenetic_cooccurrence%>% as.numeric(),
         method = "sp")

"No"

# homology



cor.test(String.ineraction.2$DELTA.score %>% as.numeric(),
         String.ineraction.2$homology %>% as.numeric(),
         method = "sp")

"No"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$homology%>% as.numeric(),
         method = "sp")

"No"

#coexpression

data.2 <- String.ineraction.2 %>% filter(coexpression>0)


# ggplot(String.ineraction.2)+
#   geom_point(aes(x=DELTA.score,y=experimentally_determined_interaction))


cor.test(data.2$DELTA.score%>% as.numeric(),
         data.2$coexpression %>% as.numeric(),
         method = "sp")

"Yes"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$coexpression%>% as.numeric(),
         method = "sp")
"NO"

#"experimentally_determined_interaction"
data.2 <- String.ineraction.2 %>% filter(experimentally_determined_interaction>0)

cor.test(data.2$DELTA.score %>% as.numeric(),
         data.2$experimentally_determined_interaction %>% as.numeric(),
         method = "sp")

"Yes"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$experimentally_determined_interaction %>% as.numeric(),
         method = "sp")

"No"

# database_annotated
data.2 <- String.ineraction.2 %>% filter(database_annotated>0)


cor.test(String.ineraction.2$DELTA.score %>% as.numeric(),
         String.ineraction.2$database_annotated %>% as.numeric(),
         method = "sp")

"Yes"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$database_annotated %>% as.numeric(),
         method = "sp")

"Yes"



# automated_textmining
data.2 <- String.ineraction.2 %>% filter(database_annotated>0)
cor.test(String.ineraction.2$DELTA.score %>% as.numeric(),
         String.ineraction.2$automated_textmining %>% as.numeric(),
         method = "sp")

"Yes"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$automated_textmining %>% as.numeric(),
         method = "sp")

"No"


# combined_score

cor.test(String.ineraction.2$DELTA.score %>% as.numeric(),
         String.ineraction.2$combined_score %>% as.numeric(),
         method = "sp")

"Yes"

cor.test(-log10(String.ineraction.3$neg.log10.pvalue %>% as.numeric()),
         String.ineraction.3$combined_score %>% as.numeric(),
         method = "sp")

"Yes"
# 
# 
# String.ineraction.2.summarise <- 
#   String.ineraction.2 %>% 
#   group_by(Sub.group) %>% 
#   summarise(Score.mean =mean(as.numeric(DELTA.score)),
#             coexpress.mean= mean(coexpression),
#             coexpress.se = std(coexpression),
#             EDI.mean = mean(experimentally_determined_interaction),
#             EDI.se = std(experimentally_determined_interaction),
#             DA.mean = mean(database_annotated),
#             DA.se = std(database_annotated),
#             AT.mean = mean(automated_textmining),
#             AT.se = std(automated_textmining),
#             CS.mean = mean(combined_score),
#             CS.se = std(combined_score)
#   )
#             
# String.ineraction.3.summarise <- 
#   String.ineraction.3 %>% 
#   group_by(Sub.group) %>% 
#   summarise(Score.mean = mean(as.numeric(DELTA.score)),
#             coexpress.mean= mean(coexpression),
#             coexpress.se = std(coexpression),
#             EDI.mean = mean(experimentally_determined_interaction),
#             EDI.se = std(experimentally_determined_interaction),
#             DA.mean = mean(database_annotated),
#             DA.se = std(database_annotated),
#             AT.mean = mean(automated_textmining),
#             AT.se = std(automated_textmining),
#             CS.mean = mean(combined_score),
#             CS.se = std(combined_score)
#   )            
# 



label_size =8

###

  
S.data <- 
  String.ineraction  %>%  
  #filter(as.numeric(experimentally_determined_interaction)>0) %>% 
  dplyr::arrange(as.numeric(DELTA.score)) %>% 
  #mutate(Index=1:nrow(.)) %>% 
  #filter(Index %in% sample(1:nrow(.),nrow(.)-nrow(.)%%10)) %>% 
  dplyr::mutate(Index=1:nrow(.)) 

S.data$Sub.group <- sapply(S.data$Index, function(x){
  
  min(which(x <= quantile(1:nrow(S.data),c(1:10)/10)))
  
})
  
S.data.group <- 
  S.data %>% 
  group_by(Sub.group)  %>% 
  summarise(Score.mean = mean(as.numeric(DELTA.score)),
           "EDI.mean" = mean(as.numeric(experimentally_determined_interaction)),
           "EDI.se" = std(as.numeric(experimentally_determined_interaction))
           )

Test <- cor.test(S.data$experimentally_determined_interaction,S.data$DELTA.score %>% as.numeric(),method="spearman")

label <- substitute(paste(rho, " = ", estimate,","," P = ", pvalue),
                    list(estimate = signif(Test$estimate, 2), pvalue = signif(Test$p.value, 2)))

p.Str <-
  ggplot(data=S.data.group)+
  geom_errorbar(aes(x=Score.mean,
                    ymin=EDI.mean-EDI.se,
                    ymax=EDI.mean+EDI.se)
  )+
  #geom_point(data=S.data,aes(x=as.numeric(DELTA.score),y=experimentally_determined_interaction,group=Sub.group,color=Sub.group))+
  geom_point(aes(x=Score.mean,y=EDI.mean,group=Sub.group))+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("Experimentally determined interaction")+
  xlab("DELTA score")

p.Str.1 <- ggdraw(p.Str) + draw_label(label, .5,0.9,size=label_size) 
p.Str.1



###


S.data <- 
  String.ineraction  %>%  
  filter(as.numeric(coexpression)>0) %>% 
  filter(!is.na(DELTA.score)) %>% 
  arrange(as.numeric(DELTA.score)) %>% 
  #mutate(Index=1:nrow(.)) %>% 
  #filter(Index %in% sample(1:nrow(.),nrow(.)-nrow(.)%%10)) %>% 
  mutate(Index=1:nrow(.)) 

S.data$Sub.group <- sapply(S.data$Index, function(x){
  
  min(which(x <= quantile(1:nrow(S.data),c(1:10)/10)))
  
})

S.data.group <- 
  S.data %>% 
  group_by(Sub.group)  %>% 
  summarise(Score.mean = mean(as.numeric(DELTA.score)),
            "coexpression.mean" = mean(as.numeric(coexpression)),
            "coexpression.se" = std(as.numeric(coexpression))
  )

Test <- cor.test(S.data$coexpression,S.data$DELTA.score %>% as.numeric(),method="spearman")

label <- substitute(paste(rho, " = ", estimate,","," P = ", pvalue),
                    list(estimate = signif(Test$estimate, 2), pvalue = signif(Test$p.value, 2)))

p.Str <-
  ggplot(data=S.data.group)+
  geom_errorbar(aes(x=Score.mean,
                    ymin=coexpression.mean-coexpression.se,
                    ymax=coexpression.mean+coexpression.se)
  )+
  #geom_point(data=S.data,aes(x=as.numeric(DELTA.score),y=coexpression,group=Sub.group,color=Sub.group))+
  geom_point(aes(x=Score.mean,y=coexpression.mean,group=Sub.group))+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("Coexpression")+
  xlab("DELTA score")

p.Str.2 <- ggdraw(p.Str) + draw_label(label, .7, .15,size=9) 
p.Str.2



###


S.data <- 
  String.ineraction.2 %>%  
  filter(as.numeric(database_annotated)>0) %>% 
  arrange(as.numeric(DELTA.score)) %>% 
  #mutate(Index=1:nrow(.)) %>% 
  #filter(Index %in% sample(1:nrow(.),nrow(.)-nrow(.)%%10)) %>% 
  mutate(Index=1:nrow(.))

S.data$Sub.group <- sapply(S.data$Index, function(x){
  
  min(which(x <= quantile(1:nrow(S.data),c(1:10)/10)))
  
})

S.data.group <- 
  S.data %>% 
  group_by(Sub.group)  %>% 
  summarise(Score.mean = mean(as.numeric(DELTA.score)),
            "DA.mean" = mean(as.numeric(database_annotated)),
            "DA.se" = std(as.numeric(database_annotated))
  )

Test <- cor.test(S.data$database_annotated,S.data$DELTA.score %>% as.numeric(),method="spearman")

label <- substitute(paste("Spearman ", rho, " = ", estimate,","," P = ", pvalue),
                    list(estimate = signif(Test$estimate, 2), pvalue = signif(Test$p.value, 2)))

p.Str <-
  ggplot(data=S.data.group)+
  geom_errorbar(aes(x=Score.mean,
                    ymin=DA.mean-DA.se,
                    ymax=DA.mean+DA.se)
  )+
  #geom_point(data=S.data,aes(x=as.numeric(DELTA.score),y=coexpression,group=Sub.group,color=Sub.group))+
  geom_point(aes(x=Score.mean,y=DA.mean,group=Sub.group))+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("Database annotated")+
  xlab("DELTA score")

p.Str.3<- ggdraw(p.Str) + draw_label(label, .8, .9,size=label_size) 
p.Str.3

###


S.data <- 
  String.ineraction.2 %>%  
  filter(as.numeric(automated_textmining)>0) %>% 
  arrange(as.numeric(DELTA.score)) %>% 
  #mutate(Index=1:nrow(.)) %>% 
  #filter(Index %in% sample(1:nrow(.),nrow(.)-nrow(.)%%10)) %>% 
  mutate(Index=1:nrow(.)) 

S.data$Sub.group <- sapply(S.data$Index, function(x){
  
  min(which(x <= quantile(1:nrow(S.data),c(1:10)/10)))
  
})


S.data.group <- 
  S.data %>% 
  group_by(Sub.group)  %>% 
  summarise(Score.mean = mean(as.numeric(DELTA.score)),
            "AT.mean" = mean(as.numeric(automated_textmining)),
            "AT.se" = std(as.numeric(automated_textmining))
  )

Test <- cor.test(S.data$automated_textmining,S.data$DELTA.score %>% as.numeric(),method="spearman")

label <- substitute(paste("Spearman ", rho, " = ", estimate,","," P = ", pvalue),
                    list(estimate = signif(Test$estimate, 2), pvalue = signif(Test$p.value, 2)))

p.Str <-
  ggplot(data=S.data.group)+
  geom_errorbar(aes(x=Score.mean,
                    ymin=AT.mean-AT.se,
                    ymax=AT.mean+AT.se)
  )+
  #geom_point(data=S.data,aes(x=as.numeric(DELTA.score),y=coexpression,group=Sub.group,color=Sub.group))+
  geom_point(aes(x=Score.mean,y=AT.mean,group=Sub.group))+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("Automated textmining")+
  xlab("DELTA score")

p.Str.4 <- ggdraw(p.Str) + draw_label(label, .5, .9,size=label_size) 
p.Str.4

###


S.data <- 
  String.ineraction.2 %>%  
  filter(as.numeric(combined_score)>0) %>% 
  arrange(as.numeric(DELTA.score)) %>% 
  #mutate(Index=1:nrow(.)) %>% 
  #filter(Index %in% sample(1:nrow(.),nrow(.)-nrow(.)%%10)) %>% 
  mutate(Index=1:nrow(.)) 

S.data$Sub.group <- sapply(S.data$Index, function(x){
  
  min(which(x <= quantile(1:nrow(S.data),c(1:10)/10)))
  
})

S.data.group <- 
  S.data %>% 
  group_by(Sub.group)  %>% 
  summarise(Score.mean = mean(as.numeric(DELTA.score)),
            "CS.mean" = mean(as.numeric(combined_score)),
            "CS.se" = std(as.numeric(combined_score))
  )

Test <- cor.test(S.data$combined_score,S.data$DELTA.score %>% as.numeric(),method="spearman")

label <- substitute(paste("Spearman ", rho, " = ", estimate,","," P = ", pvalue),
                    list(estimate = signif(Test$estimate, 2), pvalue = signif(Test$p.value, 2)))

p.Str <-
  ggplot(data=S.data.group)+
  geom_errorbar(aes(x=Score.mean,
                    ymin=CS.mean-CS.se,
                    ymax=CS.mean+CS.se)
  )+
  #geom_point(data=S.data,aes(x=as.numeric(DELTA.score),y=coexpression,group=Sub.group,color=Sub.group))+
  geom_point(aes(x=Score.mean,y=CS.mean,group=Sub.group))+
  theme_classic()+
  theme(panel.background = element_rect(color="black",size = 1.5))+
  ylab("Combined score")+
  xlab("DELTA score")

p.Str.5 <- ggdraw(p.Str) + draw_label(label, .5, .9,size=label_size) 
p.Str.5


pp.Str <- cowplot::plot_grid(p.Str.1,p.Str.2,p.Str.3,p.Str.4,p.Str.5,nrow = 1,labels = c("A","B","C","D","E"))
pp.Str <- cowplot::plot_grid(p.Str.1,p.Str.2,nrow = 1,labels = c("A","B","C","D","E"))

pp.Str
#ggsave("./result2/gene_network_c.pdf",width = 25,height = 4)




# 
# S.data.group <- 
#   S.data %>% 
#   group_by(Sub.group)  %>% 
#   summarise(Score.mean =mean(as.numeric(DELTA.score)),
#             coexpress.mean= mean(coexpression),
#             coexpress.se = std(coexpression),
#             EDI.mean = mean(experimentally_determined_interaction),
#             EDI.se = std(experimentally_determined_interaction),
#             DA.mean = mean(database_annotated),
#             DA.se = std(database_annotated),
#             AT.mean = mean(automated_textmining),
#             AT.se = std(automated_textmining),
#             CS.mean = mean(combined_score),
#             CS.se = std(combined_score)
#   )


