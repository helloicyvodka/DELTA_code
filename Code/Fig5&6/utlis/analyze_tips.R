library(plyr)
library(dplyr)
library(pipeR)
library(ggplot2)
library(pipeR)
library(data.table)
library(parallel)




mut.df.list.tip <- mclapply(list.files("~/2017-2018/tree/tree_fate_change/data/data_from_digital/mutated_for_DELTA/mutated_for_DELTA/",full.names = T),
                            function(x)as.data.frame(fread(x)))

mut.df.list.analysis <- mut.df.list.tip %>% sapply(nrow) 

mut.df.list.analysis.300 <- mut.df.list.analysis[ mut.df.list.analysis>300]

ggplot()+geom_histogram(aes(mut.df.list.analysis.300),binwidth = 1)

ggplot()+geom_histogram(aes(mut.df.list.analysis),binwidth = 1)

wt.df.tip <- as.data.frame(fread("~/2017-2018/tree/tree_fate_change/data/data_from_digital/WT_merge/wt_for_DELTA.alm"))

ggplot()+geom_histogram(aes(mut.df.list.analysis),binwidth = 1)+geom_vline(xintercept = nrow(wt.df.tip),color="red")+scale_x_continuous(breaks=seq(0,450,50))
