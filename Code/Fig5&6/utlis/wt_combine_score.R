library(plyr)
library(dplyr)
library(data.table)
library(parallel)
library(rlist)



wt <-list()
wt$file.list <- list.files("./data/data_from_digital/WT/")
wt$file.list <- wt$file.list[wt$file.list!="readme.txt"]
wt$HYP <- wt$file.list [grepl("RW10348",wt$file.list )]
wt$NEU <- wt$file.list [grepl("RW10434",wt$file.list )]
wt$PHA <- wt$file.list [grepl("RW10425",wt$file.list )]


wt.dt <- lapply(c("HYP", "NEU", "PHA"), function(m) {
  d <- lapply(wt[[m]], function(x)
    fread(x))
  names(d) <- wt[[m]]
  d
})
names(wt.dt) <- c("HYP","NEU","PHA")

# com = common

com.cell.HYP <- Reduce(intersect,lapply(wt.dt$HYP,function(x)x$cell_name))
com.cell.NEU <- Reduce(intersect,lapply(wt.dt$NEU,function(x)x$cell_name))
com.cell.PHA <- Reduce(intersect,lapply(wt.dt$PHA,function(x)x$cell_name))

com.cell <- Reduce(intersect,list(com.cell.HYP,com.cell.NEU ,com.cell.PHA ))



for(t in c("HYP","NEU","PHA")){
  for(i in 1:10){
   names(wt.dt[[t]][[i]])[names(wt.dt[[t]][[i]])=="binary_call"] <- paste0(t,"_",i)

  }
}


com.dt <- list()

for(t in c("HYP","NEU","PHA")){

    dt <- Reduce(function(x,y){merge(as.data.frame(x),as.data.frame(y)[,c(2,4)],by=c("cell_name"))},wt.dt[[t]])

    dt <- dt[,c(-2,-3,-5)]
    
    dt <- dt[ dt$cell_name %in% com.cell,]
    
    assign(paste0("com.dt.",t),dt)
    
    com.dt[[t]] <- dt
    
    write.table(dt,file=paste0("../WT_merge/",t,".txt"),quote = F,row.names = F)

    rm(dt)
}





for(t in c("HYP","NEU","PHA")){
  
    dt <- com.dt[[t]]
  
  dt[[paste0("all.",t)]] <- sapply(1:nrow(dt),function(x){
    
    r <- dt[x,2:11]
    
    if(length(r[r==1])>length(r[r==0])){
      
      y<-1
    }else{
      
      y<-0
    }
    
    y
    
  })
  
  com.dt[[t]] <-  dt
  
  rm(dt)
  
}



com.expr <-  Reduce(function(x,y){merge(as.data.frame(x),as.data.frame(y)[c(1,12)],by=c("cell_name"))},com.dt)

com.expr <- com.expr[,c(-c(2:11))]

com.expr$num.marker <- sapply(1:nrow(com.expr),function(x){
  
  r <- com.expr[x,2:4]
  
  y <- length(r[r==1])
  
  y
  
})


Tis.3 <- c("HYP","NEU","PHA")



com.expr$Class<- sapply(1:nrow(com.expr),function(x){
  
  r <- com.expr[x,2:4]
  
  y <- Tis.3[r==1]
  
  if(length(y)>1){
    y <- paste(y,collapse = "_")
  }
  
  if(length(y)==0){
    
    y<- "NOEXPR"
      
  }
  
  y
  
})




# HYP=1, NEU=2, PHA=4
com.expr$Class.num<- sapply(1:nrow(com.expr),function(x){
  
  n <- c(1,2,4)
  r <- com.expr[x,2:4]
  
  y <- n[r==1]
  
  
  
  sum(y)
  
})



write.table(com.expr,file=paste0("../WT_merge/","HYP_NEU_PHA",".txt"),quote = F,row.names = F)


































the_all_combine_WT<-expand.grid(wt.HYP,wt.NEU,wt.PHA)

names(the_all_combine_WT)<-c("HYP","NEU","PHA")



the_1000_combine_WT_0<-mclapply(1:1000,function(x){
  y<-the_all_combine_WT[x,]
  the_list<-lapply(list(as.character(y[,1]),as.character(y[,2]),as.character(y[,3])),function(the_file)fread(paste0("./data_from_digital/WT/",the_file))[,c(2,4)])
  the_df<-Reduce(function(x, y) {
    full_join(x, y, by = "cell_name")
  }, the_list)
  colnames(the_df)[2:4]<-c("HYP","NEU","PHA")
  the_df<-data.frame(Lineage=LN_2_trueLN(the_df$"cell_name"),the_df)
},mc.cores=20)



the_1000_combine_WT_0<-mclapply(the_1000_combine_WT_0,function(df){
    df$cell_type<-lapply(1:nrow(df),function(row_num){
      my_get_celltype(df[row_num,])
    }) %>% unlist()
    
    df
    
  },mc.cores = 10)

the_combn_m<-combn(c(1:1000),2) %>% t()
the_combn_sample_1000<-the_combn_m[sample(1:nrow(the_combn_m),1000,replace = F),]

################################################




the_1000_combine_WT<-mclapply(the_1000_combine_WT_0,function(df){
  
  df<-filter(df,nchar(as.character(Lineage))<=8)
  
  df
  
},mc.cores = 10)



the_1000_combine_WT<-mclapply(the_1000_combine_WT,function(x){
  the_df<-x
  the_seq<-the_df$Lineage
  the_df$nb.branches<-lapply(the_seq,function(y)nb.branches(the_seq,y)) %>% unlist()
  the_df<-filter(the_df,nb.branches==0)
  the_df
},mc.cores = 10)


lapply(1:1000,function(x){
  the_df<-the_1000_combine_WT[[x]]
  the_df<-the_df[,c(1,2,6)]
  names(the_df)<-c("Lineage","Name","Class")
  write.table(the_df,paste0("./WT_score/WT_1000_combn/",x,'.alm'),quote=F,row.names = F)
  rm_tail_blank_line(paste0("./WT_score/WT_1000_combn/",x,'.alm'))
})




mclapply(1:1000,function(x){
  y1<-paste0("./WT_score/WT_1000_combn/",the_combn_sample_1000[x,1],".alm")
  y2<-paste0("./WT_score/WT_1000_combn/",the_combn_sample_1000[x,2],".alm")
  system(paste0("/mnt/data/home/ATPs/P/HSAv0.3/2018/bin/HSA2018",
                " ",
                y1,
                " ",
                y2,
                " 20 ",
                "./gene_perturbation_20171116/fate_mutation_cost.tsv l",
                " -outfile",
                paste0(" ./gene_perturbation_20171116/tmp_alml/",x,".alml")),
         intern = F)

  },mc.cores = 10)


y1<-paste0("./WT_score/WT_1000_combn/",the_combn_sample_1000[1,1],".alm")
y2<-paste0("./WT_score/WT_1000_combn/",the_combn_sample_1000[1,2],".alm")

system(paste0("/mnt/data/home/ATPs/P/HSAv0.3/2018/bin/HSA2018",
              " ",
              "-treeS",
              " ",
              y1,
              " ",
              "-treeT",
              " ",
              y1,
              " ",
              "-max_target 100",
              " ",
              "-cost ./gene_perturbation_20171116/fate_mutation_cost.tsv",
              " ",
              "-method l",
              " ",
              "-outfile",
              " ",
              paste0("./gene_perturbation_20171116/tmp_alml/","1",".alml")
              ),
       intern = F)





paste0(" ./gene_perturbation_20171116/tmp_alml/",x,".alml")





tmp_files<-list.files("./WT_score/WT_1000_combn/") 
tmp_files<-tmp_files[grep(".almg",tmp_files)]

the_combine_WT_scores_local<-mclapply(tmp_files,
       
       function(x){
         
         the_file<-paste0("./WT_score/WT_1000_combn/",x)
         
         ReadAlml(the_file)
         
         
       },mc.cores = 10
       ) %>% unlist()








the_combine_WT_scores_local<-the_combine_WT_scores_local %>% as.numeric()

hist(the_combine_WT_scores_local,breaks = 20)#,ylim=c(0,140),xlim=c(640,780))









tmp_list<-lapply(the_1000_combine_WT,function(y){
  
  lapply(y$Lineage,function(x){nchar(as.character(x))}) %>% unlist() %>% table() %>% data.frame()
})




tmp_df<-tmp_list %>% as.data.frame() %>% t()
tmp_df<-tmp_df[seq(2,2000,by=2),]




##################################################################
tmp_rm_almg<-list.files("./WT_score/WT_1000_combn/",full.names = T) 
tmp_rm_almg<-tmp_rm_almg[grep(".almg",tmp_rm_almg)]
system(paste0("rm -f ",paste(tmp_rm_almg,collapse = " ")),intern = F)
###################################################################

names(tmp_df2)<-c("round_5","round_7","round_8","sum_tips")

tmp_df2$sum_tips<-lapply(1:nrow(tmp_df2),function(x){
  sum(as.vector(as.numeric(as.character(unlist(tmp_df2[x,1:3])))))
}) %>% unlist()


WT_files<-list.files("./data_from_digital/WT/",full.names = T)
lapply(WT_files,function(x){
 
  the_Lineage<-fread(x)[,1]
  the_Lineage<-the_Lineage[lapply(unlist(the_Lineage),function(i){
    
    
    
  })]
   
}
)




