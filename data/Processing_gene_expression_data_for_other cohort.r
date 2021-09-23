
source('functions/functions_evolution.r')
load("./data/genes_infor_downloaded_from_ensembl97.rdata")



get_GTEX_data<-function(file_name,gene_trans){
  rownames(gene_trans)<-genes_infor$genes_main_transcript
  GTEX_data<- read.csv(file_name)
  GTEX_data=GTEX_data[GTEX_data$X0!=0 & GTEX_data$X0 %in% genes_infor$genes_main_transcript ,]
  GTEX_data=GTEX_data[!(duplicated(GTEX_data$X0) | apply(GTEX_data[,-1],1,max)=="0.0"),]
  ## remove duplicate
  tissue_tpm=as.data.frame(apply(GTEX_data[,-1],2,as.numeric))
  rownames(tissue_tpm)<-gene_trans[GTEX_data$X0,]$ensembl_gene_id
  colnames(tissue_tpm)<-c(gsub("\\.+",".",colnames(GTEX_data[-1])))
  tissue_tpm$Gene.Type=genes_infor[rownames(tissue_tpm),"Gene.Type"]
  tissue_tpm$Gene.ID=rownames(tissue_tpm)
  all_tissue_data=tissue_tpm[which(!is.na(tissue_tpm$Gene.Type)),]
  colnames(all_tissue_data)=firstup(colnames(all_tissue_data))
  all_data_df=reshape2::melt(all_tissue_data,id.vars = c("Gene.ID","Gene.Type"))
  all_data_df$value=log2(all_data_df$value+1)
  tissue=as.character(unique(all_data_df$variable))
  #all_data_df=all_data_df[all_data_df$Gene.Type!="Oncogene",]
  all_data_non_cancer=all_tissue_data[all_tissue_data$Gene.Type=="Non-Cancer",tissue]
  all_data_TSG=all_tissue_data[all_tissue_data$Gene.Type=="TSG",tissue]
  
  all_data_df$len_3UTR_log=log10(genes_infor[all_data_df$Gene.ID,]$hg38_3UTR_length)
  
  all_data_df$mirna_3UTR_count=genes_infor[all_data_df$Gene.ID,]$mirna_3UTR_count
  all_data_df$hg38_3UTR_length=genes_infor[all_data_df$Gene.ID,]$hg38_3UTR_length
  all_data_df$mirna_3UTR_count_density=log10(all_data_df$mirna_3UTR_count/all_data_df$hg38_3UTR_length)
  all_data_df$rbp_site_binding_density2=genes_infor[all_data_df$Gene.ID,]$rbp_site_binding_density2
  all_data_df$rbp_site_binding_density2_support=genes_infor[all_data_df$Gene.ID,]$rbp_site_binding_density2_support
  all_data_df$hg38_3utr_ratio_all=genes_infor[all_data_df$Gene.ID,]$hg38_3utr_ratio_all
  
  colnames(all_data_df)<-c("Gene.ID","Gene.Type","cancer","norm_exp","len_3UTR_log","mirna_3UTR_count","hg38_3UTR_length","mirna_3UTR_count_density","rbp_site_binding_density2","rbp_site_binding_density2_support","hg38_3utr_ratio_all")
  
  return(all_data_df)
}

GTEXmedian_tpm<-get_GTEX_data("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_median_transcript_tpm.csv",genes_infor)

cv_tpm<-get_GTEX_data("/media/huangdan/hardisk0/HD/HD_promoter_evolution/primary_data/data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_cv_transcript_tpm.csv",genes_infor)

GTEXmedian_tpm<-GTEXmedian_tpm[which(! GTEXmedian_tpm$cancer %in% c("Testis")), ]
cv_tpm<-cv_tpm[which(! cv_tpm$cancer %in% c("Testis")), ]
GTEXmedian_tpm$cancer<-factor(as.character(GTEXmedian_tpm$cancer))
cv_tpm$cancer<-factor(as.character(cv_tpm$cancer))


rownames(GTEXmedian_tpm)<-paste(GTEXmedian_tpm[,1],GTEXmedian_tpm[,3],sep="_")
rownames(cv_tpm)<-paste(cv_tpm[,1],cv_tpm[,3],sep="_")
GTEXmedian_tpm$norm_exp_cv<-cv_tpm[rownames(GTEXmedian_tpm),]$norm_exp



GTEXmedian_tpm_high<-GTEXmedian_tpm[! GTEXmedian_tpm$norm_exp >1, ]

exp_norm_TSMI_RNA<-aggregate(GTEXmedian_tpm$norm_exp,list(GTEXmedian_tpm$Gene.ID),get_TSEI_log2)
rownames(exp_norm_TSMI_RNA)<-exp_norm_TSMI_RNA[,1]
genes_infor$TSMI_tissue_specificity_RNA=exp_norm_TSMI_RNA[genes_infor$ensembl_gene_id,2]


print("validate_microarray data")

load("/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE53000.rdata")

genes_microarray_kidney=genes_infor_exp
genes_microarray_kidney$cancer=rep("ccRCC_GSE53000_array",length(genes_microarray_kidney$ensembl_gene_id))
load("/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE103512.rdata")

genes_microarray_4_cancer=microarry_data
genes_microarray_4_cancer$cancer=paste(genes_microarray_4_cancer$cancer,"_GSE103512_array",sep="")
genes_microarray_4_cancer$cancer=gsub("cancer.type..","",genes_microarray_4_cancer$cancer)


load("/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE86958.rdata")
genes_rna_exp_lung=genes_infor_exp_rna
genes_rna_exp_lung$cancer=rep("Lung_IMA_GSE86958_RNA",length(genes_rna_exp_lung$ensembl_gene_id))
load("/media/huangdan/hardisk0/HD/HD_promoter_evolution/result/data/genes_infor_exp_GSE137327.rdata")
genes_rna_exp_colon=genes_infor_exp_rna
genes_rna_exp_colon$cancer=rep("CRC_GSE137327_RNA",length(genes_rna_exp_colon$ensembl_gene_id))


names_temp=colnames(genes_microarray_4_cancer)
all_data_other_cohort<-rbind(genes_microarray_kidney[,names_temp],
                             genes_microarray_4_cancer[,names_temp],
                             genes_rna_exp_colon[,names_temp],
                             genes_rna_exp_lung[,names_temp])

all_data_other_cohort=all_data_other_cohort[which(!is.na(all_data_other_cohort$norm_exp)),] 






save(GTEXmedian_tpm,all_data_other_cohort,file="./data/processing_validation_datasets.RData")

