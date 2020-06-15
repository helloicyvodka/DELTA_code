
# 1. import packages
library(ggplot2)
library(ggtree)
library(ggvita)


create.alm <- function(x,file.name){
  
  dt <- 
    data.frame(Lineage=names(x)) %>% 
    mutate(Name=Lineage,
           Class=sapply(Lineage,function(i)x[as.character(i)]))
  
  write.alm(dt,file.name = file.name)
  
}

# 2. testing 

create.alm(c(
  "000" = "a",
  "001" = "b",
  "010" = "c",
  "011" = "a",
  "100" = "b",
  "101" = "c",
  "110" = "a",
  "111" = "b"
),"treeS.alm")

create.alm(c(
  "00" = "b",
  "01" = "c",
  "10" = "a",
  "11" = "c"
),"treeQ.alm")


#3

create.alm(c(
  "000" = "a",
  "001" = "b",
  "01" = "c",
  "1" = "a"
),"treeQ.alm")

create.alm(c(
  "0000" = "b",
  "0001" = "a",
  "001" = "c",
  "010" = "c",
  "011" = "a",
  "1" = "b"
),"treeS.alm")



DELTA(
  treeS="./treeQ.alm",
  treeT="./treeS.alm",
  cost ="./cost.tsv",
  all = F,
  prune = 1,
  outfile ="QS.alml",
  method = "l",
  max_target = 2
)

DELTA(
  treeS="./treeQ.alm",
  treeT="./treeS.alm",
  cost ="./cost.tsv",
  all = F,
  prune = 1,
  outfile ="QS.almg",
  method = "g"
)

E1 <- readal(
  fileS="./treeQ.alm",
  fileT="./treeS.alm",
  cost ="./cost.tsv",
  all = F,
  prune = 1,
  outfile ="QS.alml",
  method = "l"
  
)

pdf("./result.pdf",paper = "a4",onefile = T)
parallel::mclapply(1:2,function(x){
  
  ggvita(E1,x)%++%
    geom_tippoint(aes(fill=I(tip.fill)),shape=21,size=8,color="NA")%++%
    geom_tiplab(aes(label=node.seq.ori),hjust=.5,color="white")
  
})
dev.off()



E2 <- readal(
  fileS="./treeQ.alm",
  fileT="./treeS.alm",
  cost ="./cost.tsv",
  all = F,
  prune = 1,
  outfile ="QS.almg",
  method = "g"
  
)

ggvita(E2,1)%++%
  geom_tippoint(aes(fill=I(tip.fill)),shape=21,size=8,color="NA")%++%
  geom_tiplab(aes(label=node.seq.ori),hjust=.5,color="white")
