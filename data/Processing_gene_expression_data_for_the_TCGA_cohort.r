#############################3
#The processed gene expression data is stored in the get_exp_data_1000_100_indepentant.rdata


options(stringsAsFactors = FALSE)
require(plyr)
library(ggpubr)
library(ggsignif)
library(cocor)
source('functions/functions_evolution.r')

get_gene_cv<-function(data_x,list_x){
  temp_data<-aggregate(data_x,list_x,function(x){mean(x,na.rm=T)})
  temp_aver<-apply(temp_data[,-1],1,get_cv)
  names(temp_aver)<-temp_data[,1]
  return(temp_aver)
}

get_gene_cv_log<-function(data_x,list_x){
  temp_data<-aggregate(data_x,list_x,function(x){mean(x,na.rm=T)})
  temp_aver<-apply(temp_data[,-1],1,get_cv_log)
  names(temp_aver)<-temp_data[,1]
  return(temp_aver)
}

get_cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
get_cv_log<-function(x){sqrt(exp(sd(x,na.rm =T)^2)-1)}
get_data1<-function(file,n){
  cancer_in <- read.csv(file, stringsAsFactors=FALSE)
  cancer_in_data=cancer_in[,-c(1:n)]
  cancer_in_infor=cancer_in[,c(1:n)]
  # if n=4 , this data is expression data then the value is log2 ;because many data's value is 0 then all data add 1
  if(n==4){cancer_in_data=log2(cancer_in_data+1)}
  rownames(cancer_in_data)<-cancer_in$transcript_id
  rownames(cancer_in_infor)<-cancer_in$transcript_id
  sample_infor<-colnames(cancer_in_data)
  sample_id=ldply(strsplit(sample_infor,".",fixed=TRUE))
  print(table(sample_id[,4]))
  normal_cancer=cancer_in_data[,sample_id[,4] %in% c("11A",'11B')]
  colnames(normal_cancer)<-substr(colnames(normal_cancer),1,12)
  tumor_cancer=cancer_in_data[,sample_id[,4] %in% c("01A",'01B')]
  colnames(tumor_cancer)<-substr(colnames(tumor_cancer),1,12)
  #tumor_cancer=get_median_sample(tumor_cancer) # take median value for duplicate samples
  #normal_cancer=get_median_sample(normal_cancer)
  # jUST RANDOMLY SELECT THE FIRST ONE 
  print("some patients have multiple samples, just select the first sample")
  print("tumor")
  print(table(duplicated(colnames(tumor_cancer))))
  print("norm")
  print(table(duplicated(colnames(normal_cancer))))
  overlap_sample=unique(colnames(normal_cancer)[colnames(normal_cancer) %in% colnames(tumor_cancer)])
  print("Number of overlapped samples ")
  print(length(overlap_sample))
  tumor_cancer=tumor_cancer[,overlap_sample]
  normal_cancer=normal_cancer[,overlap_sample]
  diff=log2(tumor_cancer/normal_cancer)
  if(n==4){diff=tumor_cancer-normal_cancer} # expression data
  gene_norm_aver<-apply(normal_cancer,1,function(x){mean(x,na.rm=T)})
  names(gene_norm_aver)<-rownames(normal_cancer)
  gene_tumor_aver<-apply(tumor_cancer,1,function(x){mean(x,na.rm=T)})
  names(gene_tumor_aver)<-rownames(tumor_cancer)
  diff_pvalue<-get_diff_sig(normal_cancer,tumor_cancer)
  all_result=list(diff=diff,norm=normal_cancer,tumor=tumor_cancer,norm_aver=gene_norm_aver,tumor_aver=gene_tumor_aver[names(gene_norm_aver)],diff_pvalue=diff_pvalue,gene_infor=cancer_in_infor)
  return(all_result)
}

get_diff_sig<-function(data_n,data_t){
  len_x=length(data_n[,1])
  p_value=c()
  for(i in 1:len_x){
    x <- try(wilcox.test(as.numeric(data_n[i,]),as.numeric(data_t[i,]),paired = T),TRUE)
    if (class(x) == "try-error") {
      p_value=c(p_value,NA)
    }else{
      p_value=c(p_value,x$p.value)
    }
  }
  names(p_value)<-rownames(data_n)
  return(p_value)
}

get_gene_aver<-function(data_x,list_x){
  temp_data<-aggregate(data_x,list_x,function(x){mean(x,na.rm=T)})
  temp_aver<-apply(temp_data[,-1],1,function(x){mean(x,na.rm=T)})
  names(temp_aver)<-temp_data[,1]
  return(temp_aver)
}

get_gene_aver2<-function(data_x,list_x){
  temp_data<-aggregate(data_x,list_x,function(x){mean(x,na.rm=T)})
  temp_aver<-temp_data[,-1]
  names(temp_aver)<-temp_data[,1]
  return(temp_aver)
}



get_all_data_independent<-function(cancers,genes_infor){
  all_result=c()
  sample_list=c()
  for(cancer in cancers){
    exp_data_temp=get_data1(paste("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/expression_data/",cancer,"_merged_select_transcript.txt",sep=""),4)
    exp_samples=colnames(exp_data_temp$diff)
    exp_trans<-names(exp_data_temp$norm_aver)[exp_data_temp$norm_aver>0 | exp_data_temp$tumor_aver>0]
    print(table(exp_trans %in%  genes_infor$genes_main_transcript))
    diff_exp<-get_gene_aver(exp_data_temp$diff,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    norm_exp<-get_gene_aver(exp_data_temp$norm,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    tumor_exp<-get_gene_aver(exp_data_temp$tumor,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    diff_exp_cv<-get_gene_cv(exp_data_temp$diff,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    norm_exp_cv<-get_gene_cv_log(exp_data_temp$norm,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    tumor_exp_cv<-get_gene_cv_log(exp_data_temp$tumor,list(exp_data_temp$gene_infor[,"ensembl_gene_id"]))
    
    
    all_gene_union=unique(exp_data_temp$gene_infor[,"ensembl_gene_id"])
    all_data=cbind(diff_exp[all_gene_union],norm_exp[all_gene_union],tumor_exp[all_gene_union], diff_exp_cv[all_gene_union],norm_exp_cv[all_gene_union],tumor_exp_cv[all_gene_union])
    rownames(all_data)<-all_gene_union
    colnames(all_data)<-c("diff_exp","norm_exp","tumor_exp", "diff_exp_cv","norm_exp_cv","tumor_exp_cv")
    all_result<-rbind(all_result,cbind(cancer,all_data))
    print(head(cbind(cancer,all_data)))
    sample_list<-rbind(sample_list,cbind(cancer,"exp",exp_samples))
  }
  return(list(all_result=all_result,sample_list=sample_list))
}






load("./data/genes_infor_downloaded_from_ensembl97.rdata.rdata")


cancers=c('BLCA',"BRCA","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")
result_independent=get_all_data_independent(cancers,genes_infor )#



de_factor<-function(data_x){
  data_x<-as.data.frame(data_x)
  f_column<-is.na(as.numeric(as.character(data_x[1,])))# check the first row to see if this coloum is numeric or character
  if(length(f_column[!f_column])>1){
    data_1<-apply(data_x[,!f_column],2,function(x) as.numeric(as.character(x)))
  }else{data_1<- as.numeric(as.character(data_x[,!f_column]))}
  
  if(length(f_column[f_column])>1){
    data_2<-apply(data_x[,f_column],2,function(x) as.character(x))
  }else{data_2<- as.character(data_x[,f_column])}
  
  data_x[,f_column]<-data_2
  data_x[,!f_column]<-data_1
  return(data_x)
}

sample_num_all=result_independent$sample_list

result_all_indepent=data.frame(de_factor(result_independent$all_result),genes_infor[substr(rownames(result_independent$all_result),1,15),])
save(result_all_indepent,sample_num_all,file="./data/get_exp_data_1000_100_indepentant.rdata")



