#########
# 根据http://www.digital-development.org/ 上的图确定13 foundercells 和他们的末端 tipcell

founder_cell_13<-c("ABala","ABalp","ABara","ABarp","ABpla","ABplp","ABpra","ABprp","MS","E","C","D","P4")
founder_cell_13<-data.frame(cell_Lineage=founder_cell_13,tip_cell_level=c(rep(5,8),6,4,5,3,3))
founder_cell_13<-mutate(founder_cell_13,tip_cell_num=2^(tip_cell_level))
founder_cell_13<-data.frame(Lineage=LN_2_trueLN(as.character(founder_cell_13$cell_Lineage)) %>% unlist(),founder_cell_13)

combine_after_founder<-function(x){
  expand.grid(rep(list(0:1),x)) %>% apply(.,1,function(x)paste0(x,collapse = ""))
}


# 所有的tip cell 的组合
all_tip_cells<-lapply(1:nrow(founder_cell_13),function(x){
  x<-founder_cell_13[x,]
  tmp_all<-combine_after_founder(x$tip_cell_level)
  paste(x$Lineage,tmp_all,sep = "")
})

# total 384 tip cells at the bottom of the tree
length(all_tip_cells %>% unlist())
#384
list.mapv(all_tip_cells,length(.))
##result: [1] 32 32 32 32 32 32 32 32 64 16 32  8  8
## correct!


########
# subset the pertubation dfs

gene_perturbation_3t_tip_cell<-lapply(gene_perturbation_3t, function(x){
  subset(x,Lineage %in% unlist(all_tip_cells))
})

list.mapv(gene_perturbation_3t_tip_cell,nrow(.)) %>% table()
# 0  10  14  18  20  28  30  34  36  42  44  46  52  54  56  58  62  64  66  70  72  74  78  80  84 
# 1   1   2   1   1   1   1   1   1   1   3   1   2   1   1   1   1   1   1   1   2   1   1   1   4 
# 88  90  92  96 106 112 116 124 126 128 130 134 136 138 140 142 144 146 148 150 152 154 156 162 164 
# 1   2   2   1   3   2   1   2   1   2   2   2   1   3   1   3   2   1   1   1   1   3   2   2   1 
# 166 168 170 172 174 176 178 180 182 184 186 188 192 194 196 198 200 204 206 208 210 212 214 216 218 
# 2   2   1   1   3   3   1   4   2   1   1   1   2   2   3   2   4   2   6   4   5   4   5   2   2 
# 222 224 226 228 230 232 236 238 240 242 244 246 248 250 256 260 274 276 292 312 314 322 330 338 350 
# 2   3   4   2   1   2   1   2   1   1   1   2   2   2   1   2   1   1   1   2   1   1   1   2   1 
# 352 356 360 362 364 366 368 370 372 374 376 378 380 382 384 
# 1   3   3   2   2   1   2   1   1   1   1   1   2   1   1 


# try tip_cell_level = 2
all_tip_cells_2<-lapply(1:nrow(founder_cell_13),function(x){
  x<-founder_cell_13[x,]
  tmp_all<-combine_after_founder(3)
  paste(x$Lineage,tmp_all,sep = "")
})


gene_perturbation_3t_tip_cell<-lapply(gene_perturbation_3t, function(x){
  subset(x,nchar(Lineage)==10)
})

list.mapv(gene_perturbation_3t_tip_cell,nrow(.)) %>% table()





# tipcell_num  48  50  52 
# count        134  13  57 



system(paste0("/mnt/data/home/ATPs/P/HSAv0.3/2018/bin/HSA2018"," ","./gene_perturbation_20171116/tmp_alml/copy_wt.alm"," ","./gene_perturbation_20171116/tmp_alml/copy_wt.alm"," ","./gene_perturbation_20171116/fate_mutation_cost.tsv g"),intern = T)





