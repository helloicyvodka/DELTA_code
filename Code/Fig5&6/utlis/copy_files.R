


lapply(list.files("./data_from_digital/WT",full.names=T)[1:30],
       
       function(file){
         df<-fread(file,header=T)
         df<-data.frame(Lineage=LN_2_trueLN(df$cell_name),df)
         write.table(df,file,quote=F,row.names=F)
       }
       
       
       
       )


the_file<-list.files("./data_from_digital/Cell_Lineage_pertution")
the_file_1<-the_file[grepl("pal-1",the_file,ignore.case = T)]

the_file_2<-the_file[grepl("hda-1",the_file,ignore.case = T)]

the_filesss<-c(the_file_1,the_file_2)


lapply(the_filesss,
       
       function(f){
         df<-fread(paste0("/mnt/data/home/yuanmeng/2017-2018/tree/tree_fate_change/data_from_digital/Cell_Lineage_pertution/",f),header=T)
         df<-data.frame(Lineage=LN_2_trueLN(df$cell_name),df)
         write.table(df,paste0("./data_from_digital/pal-1_hda-1/",f),quote=F,row.names=F)
       }
)


system("mkdir ./data_from_digital/pal-1_hda-1")
system(paste0("cp ",paste(the_filesss,collapse = " ")," ./data_from_digital/pal-1_hda-1"))
system("zip -r digital_data.zip ./data_from_digital/pal-1_hda-1 ./data_from_digital/WT")
