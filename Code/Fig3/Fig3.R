
# import packages
require(dplyr)
require(ggplot2)
require(data.table)
setwd("~/2017-2018/tree/C.elegans_VS_C.elegans/data")


# cost files and hyper-parameters
cost <- 
  c("./matrix_fun.csv",
    "./cost.common.m8u8.tsv",
    "./cost.tsv",
    "/mnt/data/home/ATPs/W/Celegans/Lineages/costMatch10otherMinus100.tsv",
    "~/2017-2018/tree/ReDrawResult/cost-same-m10u2.tsv")[5]
prune = 1
match_loci = 10
mismatch_loci = 1
test = 100
outfile = "/mnt/data/home/yuanmeng/2017-2018/tree/Fig.4B/20181128-max10000-test1000-same-m10u2pr1.txt"


## Symmetric pairs score, pvalues, and match lengths.

Sym.pair <- read.csv("~/2017-2018/tree/C.elegans_VS_C.elegans/data/symmetry-pairs.csv",stringsAsFactors = F,colClasses = "character")

Sym.pair.2 <- 
  union(
    Sym.pair[,c(1,2)],
    data.frame(Sister1=Sym.pair$Sister2,Sister2=Sym.pair$Sister1,stringsAsFactors = F)
  ) %>%
  arrange(Sister1)

Sym.pair.2$Sister1_LN <- Sym.pair.2$Sister1 %>% ggvita::LN_to_Bin()
Sym.pair.2$Sister2_LN <- Sym.pair.2$Sister2 %>% ggvita::LN_to_Bin()

Sym.pair.2 <- Sym.pair.2 %>% filter(Sister1_LN>Sister2_LN)


## extract the subtrees of symmetric pairs and calculate DELTA scores and p-values of each pair

fun.alm <-
  read.table(
    "/mnt/data/home/ATPs/W/Celegans/Lineages/fun.alm",
    stringsAsFactors = F,
    header = T,
    colClasses = "character"
  )



if (T) {
  for (i in 1:nrow(Sym.pair.2)) {
    cat(i, "\n")
    r <- Sym.pair.2[i,]
    treeS <- fun.alm[fun.alm$Lineage %>% startsWith(r$Sister1_LN),]
    treeT <- fun.alm[fun.alm$Lineage %>% startsWith(r$Sister2_LN),]
    dir.create(paste0(
      "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2"
    ))
    dir.create(
      paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i
      )
    )
    
    if (length(treeS$Lineage) > 1) {
      treeS$Lineage <- treeS$Lineage %>% sub(r$Sister1_LN, "", .)
    } else{
      treeS$Lineage <- "Root"
    }
    if (length(treeT$Lineage) > 1) {
      treeT$Lineage <- treeT$Lineage %>% sub(r$Sister2_LN, "", .)
    } else{
      treeT$Lineage <- "Root"
    }
    
    ggvita::write.alm(
      treeS,
      paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i,
        "/treeS.alm"
      )
    )
    ggvita::write.alm(
      treeT,
      paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i,
        "/treeT.alm"
      )
    )
    
    ggvita::DELTA(
      treeS  = paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i,
        "/treeS.alm"
      ),
      treeT  = paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i,
        "/treeT.alm"
      ),
      cost   = cost,
      method = "g",
      outfile =  paste0(
        "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
        i,
        "/result.almg"
      ),
      prune = prune,
      test = test,
      DELTA.address = NULL
    )
  }
}


# Store the data as list

Sym.pair.list.2 <- list()

for (i in 1:nrow(Sym.pair.2)) {
  if (length(list.files(
    paste0(
      "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
      i
    )
  )) == 3) {
    Sym.pair.list.2[[i]] <-
      ggvita::readal(
        fileS  = paste0(
          "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
          i,
          "/treeS.alm"
        ),
        fileT  = paste0(
          "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
          i,
          "/treeT.alm"
        ),
        cost   = cost,
        method = "g",
        test = test,
        outfile =  paste0(
          "~/2017-2018/tree/C.elegans_VS_C.elegans/data/Symmetrical_pairs2/",
          i,
          "/result.almg"
        ),
        prune = prune
      )
  } else{
    Sym.pair.list.2[[i]] <- NA
  }
  
}


Sym.pair.2$Score <- sapply(1:31, function(i) {
  if (!is.na(Sym.pair.list.2[[i]])) {
    Sym.pair.list.2[[i]]$Score
  } else{
    NA
  }
})

Sym.pair.2$Score <- Sym.pair.2$Score %>% as.integer()

Sym.pair.2$pvalue <- sapply(1:31, function(i) {
  if (!is.na(Sym.pair.list.2[[i]])) {
    Sym.pair.list.2[[i]]$PValue$pvalue
  } else{
    NA
  }
})

Sym.pair.2$pvalue <- Sym.pair.2$pvalue %>% as.numeric()

Sym.pair.2$Match.terminal <- sapply(1:31, function(i) {
  if (!is.na(Sym.pair.list.2[[i]])) {
    length(ggvita::Find.tips(Sym.pair.list.2[[i]]$MatchS))
  } else{
    NA
  }
})

Sym.pair.2$Match.all <- sapply(1:31, function(i) {
  if (!is.na(Sym.pair.list.2[[i]])) {
    length(Sym.pair.list.2[[i]]$MatchS)+1
  } else{
    NA
  }
})


## calculate the number of prune nodes, the number of matching terminal nodes

Find_Prune_Num <- function(x){
  x <- as.character(x)
  startsWith(fun.alm$Lineage,x) %>% sum()
  
}



Find.terminals <- function (allLin) 
{
  the.allMother <- allLin %>% substr(1, nchar(.) - 1) %>% 
    unique()
  the.tips <- dplyr::setdiff(allLin, the.allMother)
  return(the.tips)
}


Sym.pair.2$NumPruneS <-
  sapply(1:31,
         function(i) {
           if(!all(is.na(Sym.pair.list.2[[i]]))){
             if(!all(is.na(Sym.pair.list.2[[i]]$PruneS))){
               ori.subtr <- paste0(Sym.pair.2[i,]$Sister1_LN %>% as.character(),
                                   Sym.pair.list.2[[i]]$PruneS)
               if(length(ori.subtr)>1){ori.subtr <- Find.terminals(ori.subtr)}
               Find_Prune_Num(ori.subtr) %>% sum()
             }else{0}
             
           }else{0} 
         })


Sym.pair.2$NumPruneT <-
  sapply(1:31,
         function(i) {
           if(!all(is.na(Sym.pair.list.2[[i]]))){
             if(!all(is.na(Sym.pair.list.2[[i]]$PruneT))){
               ori.subtr <- paste0(Sym.pair.2[i,]$Sister2_LN %>% as.character(),
                                   Sym.pair.list.2[[i]]$PruneT)
               if(length(ori.subtr)>1){ori.subtr <- Find.terminals(ori.subtr)}
               Find_Prune_Num(ori.subtr) %>% sum()
             }else{0}
             
           }else{0} 
         })




## test
# for(i in 1:31){
#   print(i)
#   if (!all(is.na(Sym.pair.list.2[[i]]))){
#              if (!all(is.na(Sym.pair.list.2[[i]]$PruneS))){
#                ori.subtr <- paste0(Sym.pair.2[i,]$Sister1_LN %>% as.character(),
#                                    Sym.pair.list.2[[i]]$PruneS)
#                if(length(ori.subtr)>1){ori.subtr <- Find.terminals(ori.subtr)}
#                Find_Prune_Num(ori.subtr) %>% sum()
#              }
#              
#            } 
#   
# }

data.table::fwrite(Sym.pair.2,file = "../SymPair2.csv")




# Fig.3B 
# pick up the wanted pairs

Sym.pair.3 <- Sym.pair.2[c(15, 20, 22, 23, 24, 25, 26, 29), ]
Sym.pair.3$Group <-
  paste0(Sym.pair.3$Sister1, "\n", Sym.pair.3$Sister2)


# calculate the first 10000 outputs and find whether the symmetric pairs are within there.

ggvita::DELTA(treeS = "/mnt/data/home/ATPs/W/Celegans/Lineages/fun.alm",
              treeT = "/mnt/data/home/ATPs/W/Celegans/Lineages/fun.alm",
              cost = cost,
              max_target = 10000,
              test = test,
              prune = prune,
              method = "l",
              all = F,
              outfile = outfile,
              DELTA.address = "/mnt/data/home/phil/acting/treeComparison/code/DELTA/bin/Release/DELTA"
)

alignments.alml <- ggvita::readal.alml(outfile)

alignments <- lapply(1:10000,function(i){
  
  
  alist <- alignments.alml[[as.character(i)]]
  
  c("id"=i,
    "Score"=alist$Score %>% as.integer(),
    "RootS"=alist$RootS %>% as.character(),
    "RootT"=alist$RootT %>% as.character(),
    "pvalue"=alist$PValue$pvalue %>% as.character(),
    "MatchLength"=length(alist$MatchS)+1,
    "log10pvalue"=(alist$PValue$pvalue %>% as.numeric()%>% log10()),
    "MatchS"=list(alist$MatchS),
    "MatchT"=list(alist$MatchT)
  )
  
}) %>% 
  Reduce(rbind,.) %>% 
  as.data.frame()

#tbl_df(alignments)


Sym.pair.4 <- Sym.pair.3
rownames(Sym.pair.4) <- 1:nrow(Sym.pair.4)
alignments.2 <- lapply(1:nrow(Sym.pair.4),function(x){
  
  r <- Sym.pair.4[x,]
  
  tb1 <- dplyr::filter(alignments,RootS==r$Sister1_LN,RootT==r$Sister2_LN)
  tb2 <- dplyr::filter(alignments,RootS==r$Sister2_LN,RootT==r$Sister1_LN)
  
  rbind(tb1,tb2)
  
})%>% 
  Reduce(rbind,.) %>% 
  as.data.frame()

#tbl_df(alignments.2)


Sym.pair.4.MatchLength <- lapply(1:nrow(Sym.pair.4),function(x){
  
  r <- Sym.pair.4[x,]
  
  #df <- alignments %>% filter(MatchLength==r$Match.Length)
  
  df <- alignments %>% 
    dplyr::filter(MatchLength>=as.numeric(r$Match.Length)*0.9 & MatchLength<=as.numeric(r$Match.Length)*1)
  
  df <- df %>% dplyr::filter(!(RootS==r$Sister1_LN & RootT==r$Sister2_LN))
  df <- df %>% dplyr::filter(!(RootS==r$Sister2_LN & RootT==r$Sister1_LN))
  #df <- df %>% filter(id>=20)
  
  
  if(nrow(df)!=0){
    data.frame(x=as.character(r$Group),
               Score=unlist(df$Score),
               MatchLength=unlist(df$MatchLength),
               RootS=unlist(df$RootS),
               RootT=unlist(df$RootT),
               stringsAsFactors = F)
  }else{
    NULL
  }
  
  
}) %>% 
  Reduce(rbind,.)

Sym.pair.4.MatchLength$Sym <- sapply(1:nrow(Sym.pair.4.MatchLength),function(x){
  
  x <- Sym.pair.4.MatchLength[x,]
  df1 <- dplyr::filter(Sym.pair.4,Sister1_LN==x$RootS,Sister2_LN==x$RootT)
  df2 <- dplyr::filter(Sym.pair.4,Sister1_LN==x$RootT,Sister2_LN==x$RootS)
  
  nrow(rbind(df1,df2))!=0
  
})

Sym.pair.4.MatchLength <- Sym.pair.4.MatchLength %>% dplyr::filter(Sym==F)
Sym.pair.4.MatchLength$ID <- 1:nrow(Sym.pair.4.MatchLength)

Sym.pair.4.MatchLength$Replicated <- sapply(1:nrow(Sym.pair.4.MatchLength),function(x){
  
  r <- Sym.pair.4.MatchLength[x,]
  tb1 <- Sym.pair.4.MatchLength %>% dplyr::filter(x==r$x,RootS==r$RootT,RootT==r$RootS)
  if(nrow(tb1)==1){
    tb1$ID > x
  }else{NA}
})


p2 <-
  Sym.pair.4.MatchLength %>% 
  dplyr::filter(Replicated==F) %>% 
  ggplot2::ggplot()+
  geom_col(data=Sym.pair.4,aes(x=Group,y=5*as.numeric(-log10(pvalue))),fill="grey",color=NA,show.legend = F)+
  geom_jitter(aes(x=x,y=as.integer(Score)),color="blue",width = 0.1,size=.01,height = 0,show.legend = F)+
  geom_point(data=Sym.pair.4,aes(x=Group,y=as.integer(Score)),color="red")+
  #geom_text(data=Sym.pair.4,aes(x=Group,y=700,label=signif(log10(as.numeric(pvalue)),digits = 3)))+
  ylab("DELTA score")+
  xlab(NULL)+
  scale_y_continuous(breaks = c(0:5)*100,sec.axis = sec_axis(~.*0.2,breaks = c(0:5)*20,name = "-log10(P)"))+
  theme_classic()
p2
ggsave(paste0("./Fig3-",strsplit(outfile,split="/",fixed = T)[[1]][9],".pdf"),plot = p2,width =20,units = "cm",height = 5)




# p3 <-
#   Sym.pair.3 %>% 
#   ggplot2::ggplot()+
#   geom_col(aes(x=Group,y=as.numeric(Score),fill=Group),color=NA,show.legend = F)+
#   ylab("Score")+
#   xlab(NULL)+
#   theme_classic()
# 
# p4 <-
#   Sym.pair.3 %>% 
#   ggplot2::ggplot()+
#   geom_col(aes(x=Group,y=as.numeric(-log10(pvalue))),fill="grey",color=NA,show.legend = F)+
#   ylab("-log10(p-value)")+
#   theme_classic()
# 
# 
# cowplot::plot_grid(p2,p4,ncol=1,rel_heights = 1)





