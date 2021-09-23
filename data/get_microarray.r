source('functions/functions_evolution.r')
load("./data/genes_infor_downloaded_from_ensembl97.rdata")

library(Biobase)
library(GEOquery)
library(limma)
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




gset <- getGEO("GSE17612", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

sample_infor=gset@phenoData@data
genes_information=gset@featureData@data
gene_data= exprs(gset)
gene_data_with_ID=gene_data[as.character(genes_information[genes_information$`Gene symbol`!="",]$ID),]
rownames(gene_data_with_ID)<-genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`
x<-do.call(rbind,strsplit(genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`,"///"))

temp_x=aggregate(gene_data_with_ID,list(as.character(x[,2])),mean)
gene_data_all=temp_x[,-1]
rownames(gene_data_all)<-temp_x[,1]

sample.tissue<-do.call(rbind,strsplit(sample_infor[colnames(gene_data_all),]$title,"_"))[,1]
#select.tissue<-table(sample.tissue)
gene_data_all.tissue<-aggregate(t(gene_data_all),list(sample.tissue),median)

#################################################

gset <- getGEO("GSE23546", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
sample_infor=gset@phenoData@data
genes_information=gset@featureData@data
gene_data= exprs(gset)
gene_data_with_ID=gene_data[as.character(genes_information[genes_information$`Gene symbol`!="",]$ID),]
rownames(gene_data_with_ID)<-genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`
x<-do.call(rbind,strsplit(genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`,"///"))
temp_x=aggregate(gene_data_with_ID,list(as.character(x[,2])),mean)
gene_data_all=temp_x[,-1]
rownames(gene_data_all)<-temp_x[,1]

cv<-apply(gene_data_all,1,get_cv_log)
med<-apply(gene_data_all,1,median)
genes_infor$cv_brain<-cv[genes_infor$hgnc_symbol]
genes_infor$med_brain<-med[genes_infor$hgnc_symbol]


save(gene_data_all,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE23546.rdata')  


##############################################3

gset <- getGEO("GSE3526", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

sample_infor=gset@phenoData@data
genes_information=gset@featureData@data
gene_data= exprs(gset)
gene_data_with_ID=gene_data[as.character(genes_information[genes_information$`Gene symbol`!="",]$ID),]
rownames(gene_data_with_ID)<-genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`
x<-do.call(rbind,strsplit(genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`,"///"))

temp_x=aggregate(gene_data_with_ID,list(as.character(x[,2])),mean)
gene_data_all=temp_x[,-1]
rownames(gene_data_all)<-temp_x[,1]

sample.tissue<-do.call(rbind,strsplit(sample_infor[colnames(gene_data_all),]$title,"_"))[,1]
#select.tissue<-table(sample.tissue)
gene_data_all.tissue<-aggregate(t(gene_data_all),list(sample.tissue),median)
gene_data_all.tissue.cv<-aggregate(t(gene_data_all),list(sample.tissue),get_cv)


genes_infor_exp_GSE3526=list(median=gene_data_all.tissue,cv=gene_data_all.tissue.cv)
save(genes_infor_exp_GSE3526,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE3526.rdata')  

##############3

gset <- getGEO("GSE53000", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

sample_infor=gset@phenoData@data
genes_information=gset@featureData@data
gene_data= exprs(gset)
gene_data_with_ID=gene_data[as.character(genes_information[genes_information$`Gene symbol`!="",]$ID),]
rownames(gene_data_with_ID)<-genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`

x<-do.call(rbind,strsplit(genes_information[genes_information$`Gene symbol`!="",]$`Gene symbol`,"///"))

temp_x=aggregate(gene_data_with_ID,list(as.character(x[,2])),mean)
gene_data_all=temp_x[,-1]
rownames(gene_data_all)<-temp_x[,1]

normal_Sample=sample_infor[sample_infor$`tissue:ch1`=="normal adult kidney",]
data_n=gene_data_all[,rownames(sample_infor[sample_infor$`tissue:ch1`=="normal adult kidney",])]
colnames(data_n)=normal_Sample$`patient identifier:ch1`
tumor_Sample=sample_infor[sample_infor$`patient identifier:ch1` %in% normal_Sample$`patient identifier:ch1` & sample_infor$`tissue:ch1` == "clear cell renal cell carcinoma",]
temp_t=gene_data_all[,rownames(tumor_Sample)]
temp_t1=aggregate(t(temp_t),list(tumor_Sample$`patient identifier:ch1`),mean)
data_t=t(temp_t1[,-1])
colnames(data_t)<-temp_t1[,1]
print(dim(data_n))
diff_sig=get_diff_sig(data_n,data_t)
N_mean=apply(data_n,1,function(x){mean(x,na.rm=T)})
N_cv=apply(data_n,1,get_cv_log)
T_mean=apply(data_t,1,function(x){mean(x,na.rm=T)})
log2_FC=T_mean-N_mean
genes=rownames(data_n)
all_result=data.frame(ens_id=genes,diff_exp_pvalue=diff_sig[genes],diff_exp=log2_FC[genes],norm_exp=N_mean[genes],tumor_exp=T_mean[genes],norm_exp_cv=N_cv[genes])
rownames(all_result)<-genes

genes_infor_exp=data.frame(genes_infor,all_result[genes_infor$hgnc_symbol,])
save(genes_infor_exp_GSE53000=genes_infor_exp,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE53000.rdata')  

gset <- getGEO("GSE103512", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))



sample_infor=gset@phenoData@data
genes_information=gset@featureData@data
gene_data= exprs(gset)
gene_data_with_ID=gene_data[genes_information[genes_information$Gene.symbol!="",]$ID,]
rownames(gene_data_with_ID)<-genes_information[genes_information$Gene.symbol!="",]$Gene.symbol
temp_x=aggregate(gene_data_with_ID,list(genes_information[genes_information$Gene.symbol!="",]$Gene.symbol),mean)
gene_data_all=temp_x[,-1]
rownames(gene_data_all)<-temp_x[,1]

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

get_data<-function(x_data,cancer,sample_infor){
  sample_infor_ids=sample_infor$geo_accession
  names(sample_infor_ids)<-sample_infor$title
  sample_normal_names=sample_infor$title[sample_infor$characteristics_ch1.1=="normal: yes" & sample_infor$characteristics_ch1==cancer]
  sample_tumor_names= gsub("_normal","",sample_normal_names)
  
  data_n=x_data[,sample_infor_ids[sample_normal_names]]
  data_t=x_data[,sample_infor_ids[sample_tumor_names]]
  diff_sig=get_diff_sig(data_n,data_t)
  N_mean=apply(data_n,1,function(x){mean(x,na.rm=T)})
  N_cv=apply(data_n,1,get_cv_log)
  T_mean=apply(data_t,1,function(x){mean(x,na.rm=T)})
  log2_FC=T_mean-N_mean
  genes=rownames(data_n)
  all_result=data.frame(cancer=rep(make.names(cancer),length(genes)),gene_symbol=genes,diff_exp_pvalue=diff_sig[genes],diff_exp=log2_FC[genes],norm_exp=N_mean[genes],tumor_exp=T_mean[genes],norm_exp_cv=N_cv[genes])
  rownames(all_result)<-genes
  return(all_result)
}


breast=get_data(gene_data_all,"cancer type: BC",sample_infor)
colon=get_data(gene_data_all,"cancer type: CRC",sample_infor)
lung=get_data(gene_data_all,"cancer type: NSCLC",sample_infor)
prostate=get_data(gene_data_all,"cancer type: PCA",sample_infor)

all_data_microarry=rbind(breast,colon,lung,prostate)

microarry_data=merge(genes_infor,all_data_microarry,by.x="hgnc_symbol",by.y="gene_symbol")
save(genes_infor_exp_GSE103512=microarry_data,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data//genes_infor_exp_GSE103512.rdata')




ensembl_NCBI_transcript_hg38 <- read.delim("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/ensembl_NCBI_transcript_hg38.txt")

es_nc=ensembl_NCBI_transcript_hg38[ensembl_NCBI_transcript_hg38$RefSeq.mRNA.ID!="",]
es_nc=es_nc[es_nc$Transcript.stable.ID %in% genes_infor$genes_main_transcript,]

#################################################################### RNA -seq 
file_x=c("GSM2311004_Case1_Tumor_fpkm.txt","GSM2311010_Case1_Normal_fpkm.txt",
         "GSM2311005_Case2_Tumor_fpkm.txt","GSM2311011_Case2_Nomal_fpkm.txt",
         "GSM2311006_Case3_Tumor_fpkm.txt","GSM2311012_Case3_Normal_fpkm.txt",
         "GSM2311007_Case4_Tumor_fpkm.txt","GSM2311013_Case4_Normal_fpkm.txt",
         "GSM2311009_Case6_Tumor_fpkm.txt","GSM2311015_Case6_Normal_fpkm.txt",
         "GSM2311008_Case5_Tumor_fpkm.txt","GSM2311014_Case5_Normal_fpkm.txt")



file_n<-substr(file_x,12,22)
all_data=genes_infor[,c("genes_main_transcript","Gene.Type")]
rownames(all_data)<-all_data[,1]
for(file_1 in file_x){
  isoforms <- read.delim(paste("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GSE86958_RAW/",file_1,sep=""))
  isoforms$tracking_id=as.character(isoforms$tracking_id)
  print(file_1)
  print(head(isoforms))
  if(substr(isoforms$tracking_id[1],1,1)=="N"){
    #     Gene.stable.ID Transcript.stable.ID RefSeq.mRNA.ID
    #58  ENSG00000176555      ENST00000319988   NM_001004725
    #78  ENSG00000141682      ENST00000316660      NM_021127
    #79  ENSG00000104903      ENST00000264824      NM_005583
    temp=ldply(strsplit(isoforms$tracking_id,".",fixed = T))
    isoforms$tracking_id=temp[,1]   
    es_nc_temp=es_nc[es_nc$RefSeq.mRNA.ID %in% isoforms$tracking_id ,]
    es_nc_temp=es_nc_temp[!duplicated(es_nc_temp$RefSeq.mRNA.ID),]
    rownames(es_nc_temp)<-es_nc_temp$RefSeq.mRNA.ID
    
    isoforms=isoforms[isoforms$tracking_id %in% es_nc$RefSeq.mRNA.ID,]
    isoforms=isoforms[!duplicated(es_nc_temp[isoforms$tracking_id,]$Transcript.stable.ID),]
    rownames(isoforms)<-es_nc_temp[isoforms$tracking_id,]$Transcript.stable.ID
  }else{
    temp=ldply(strsplit(isoforms$tracking_id,".",fixed = T))
    rownames(isoforms)=temp[,1]     
  }
  all_data=cbind(all_data,isoforms[rownames(all_data),"FPKM"])
  colnames(all_data)<-c(colnames(all_data)[1:(length(all_data[1,])-1)],substr(file_1,12,22))
}

all_data_matrix=all_data[,-c(1:2)]

types=substr(colnames(all_data_matrix),7,11)

data_n=all_data_matrix[,types %in% c("Norma","Nomal")]
data_t=all_data_matrix[,types=="Tumor"]
print(dim(data_n))
diff_sig=get_diff_sig(data_n,data_t)
N_mean=apply(data_n,1,function(x){mean(x,na.rm=T)})
N_cv=apply(data_n,1,get_cv)
T_mean=apply(data_t,1,function(x){mean(x,na.rm=T)})
log2_FC=log2(T_mean/N_mean)
genes=rownames(data_n)
all_result_rna=data.frame(ens_id=genes,diff_exp_pvalue=diff_sig[genes],diff_exp=log2_FC[genes],norm_exp=N_mean[genes],tumor_exp=T_mean[genes],norm_exp_cv=N_cv[genes])
rownames(all_result_rna)<-genes

genes_infor_exp_rna=data.frame(genes_infor,all_result_rna[genes_infor$genes_main_transcript,])

save(genes_infor_exp_GSE86958=genes_infor_exp_rna,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE86958.rdata') 



#################################################################### RNA -seq 
file_x=c("GSM4075466_C1.gene.fpkm.txt","GSM4075467_T1.gene.fpkm.txt",
         "GSM4075468_C2.gene.fpkm.txt","GSM4075469_T2.gene.fpkm.txt",
         "GSM4075470_C3.gene.fpkm.txt","GSM4075471_T3.gene.fpkm.txt",
         "GSM4075472_C4.gene.fpkm.txt","GSM4075473_T4.gene.fpkm.txt",
         "GSM4075474_C5.gene.fpkm.txt","GSM4075475_T5.gene.fpkm.txt",
         "GSM4075476_C6.gene.fpkm.txt","GSM4075477_T6.gene.fpkm.txt",
         "GSM4075478_C7.gene.fpkm.txt","GSM4075479_T7.gene.fpkm.txt",
         "GSM4075480_C8.gene.fpkm.txt","GSM4075481_T8.gene.fpkm.txt",
         "GSM4075482_C9.gene.fpkm.txt","GSM4075483_T9.gene.fpkm.txt")



file_n<-substr(file_x,12,13)
isoforms <- read.delim(paste("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GSE137327_RAW/",file_x[1],sep=""))
all_data=isoforms[,c("gene_id","Symbol")]
rownames(all_data)<-all_data[,1]
for(file_1 in file_x){
  isoforms <- read.delim(paste("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GSE137327_RAW/",file_1,sep=""))
  rownames(isoforms)=isoforms$gene_id
  print(file_1)
  all_data=cbind(all_data,isoforms[rownames(all_data),"FPKM"])
  colnames(all_data)<-c(colnames(all_data)[1:(length(all_data[1,])-1)],substr(file_1,12,13))
}

all_data_matrix_temp=all_data[!is.na(all_data$Symbol),-c(1:2)]
print(all_data_matrix_temp[1:3,1:5])
temp_data=aggregate(all_data_matrix_temp,list(all_data[!is.na(all_data$Symbol),]$Symbol),function(x){mean(x,na.rm=T)})
types=substr(colnames(all_data_matrix_temp),1,1)
all_data_matrix=temp_data[,-1]
rownames(all_data_matrix)<-temp_data[,1]
data_n=all_data_matrix[,types %in% c("C")]
data_t=all_data_matrix[,types=="T"]
print(dim(data_n))
diff_sig=get_diff_sig(data_n,data_t)
N_mean=apply(data_n,1,function(x){mean(x,na.rm=T)})
T_mean=apply(data_t,1,function(x){mean(x,na.rm=T)})
N_cv=apply(data_n,1,get_cv)
log2_FC=log2(T_mean/N_mean)
genes=rownames(data_n)
all_result_rna=data.frame(ens_id=genes,diff_exp_pvalue=diff_sig[genes],diff_exp=log2_FC[genes],norm_exp=N_mean[genes],tumor_exp=T_mean[genes],norm_exp_cv=N_cv[genes])
rownames(all_result_rna)<-genes

genes_infor_exp_rna=data.frame(genes_infor,all_result_rna[genes_infor$hgnc_symbol,])
print(genes_infor_exp_rna)
save(genes_infor_exp_137327=genes_infor_exp_rna,file='/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE137327.rdata') 

