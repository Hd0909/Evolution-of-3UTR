
source('functions/functions_evolution.r')
load("./data/processing_validation_datasets.RData")
load("./data/get_exp_data_1000_100_indepentant.rdata")

load("./data/genes_infor_downloaded_from_ensembl97.rdata")

print(head(result_all_indepent))
print(str(result_all_indepent))
all_data_tsg=result_all_indepent
gene_age_infor <- read.delim("./data/gene_age_infor.txt")
rownames(gene_age_infor)<-gene_age_infor$Age.class
genes_infor$Origin.time..million.years.ago..<-gene_age_infor[genes_infor$Age_class,"Origin.time..million.years.ago.."]
##### dividing the genes into different classes based on the gene ages
cut_x=c(0,7,14,18,23,27)
cut_lable=c()
for(i in 2:length(cut_x)){
  x=cut_x[i-1]
  if(x==0){x=1}
  cut_lable<-c(cut_lable,paste("(",gene_age_infor[x,4]," - ",gene_age_infor[cut_x[i],4],"]",sep=""))
}
cut_lable[5]=">1119.25"

genes_infor<-genes_infor %>% mutate(age_type2=cut(genes_infor$Age_class, breaks=cut_x,labels = cut_lable))
rownames(genes_infor)<-genes_infor$ensembl_gene_id
all_data_tsg$age_type2=genes_infor[all_data_tsg$ensembl_gene_id,"age_type2"]



rownames(genes_infor)<-genes_infor$ensembl_gene_id







all_data_tsg$hg38_3utr_to_cds=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_3utr_to_cds
all_data_tsg$genes_trans_length_log<-log10(genes_infor[all_data_tsg$ensembl_gene_id,]$genes_trans_length)
all_data_tsg$cds_ratio_all<-log10(genes_infor[all_data_tsg$ensembl_gene_id,]$cds_ratio_all)

all_data_tsg$hg38_3utr_ratio_non5UTR=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_3utr_ratio_non5UTR
all_data_tsg$log_ratio_3UTR_5UTR=genes_infor[all_data_tsg$ensembl_gene_id,]$log_ratio_3UTR_5UTR
all_data_tsg$mirna_binding_3utr_density=genes_infor[all_data_tsg$ensembl_gene_id,]$mirna_binding_3utr_density
all_data_tsg$mirna_3UTR_count_density=genes_infor[all_data_tsg$ensembl_gene_id,]$mirna_3UTR_count_density

#all_data_tsg$structure_length_log=genes_infor[all_data_tsg$ensembl_gene_id,]$structure_length_log
all_data_tsg$hg38_5utr_to_cds=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_5utr_to_cds
all_data_tsg$hg38_3utr_ratio_all=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_3utr_ratio_all
all_data_tsg$hg38_5utr_ratio_all=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_5utr_ratio_all
all_data_tsg$hg38_3UTR_length=genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_3UTR_length
all_data_tsg$rbp_site_count_density=genes_infor[all_data_tsg$ensembl_gene_id,]$rbp_site_count_density
all_data_tsg$len_3UTR_log_bins=genes_infor[all_data_tsg$ensembl_gene_id,]$len_3UTR_log_bins
all_data_tsg$ratio_all_3UTR_log_bins=genes_infor[all_data_tsg$ensembl_gene_id,]$ratio_all_3UTR_log_bins
all_data_tsg$hg38_3UTR_length_log=log10(genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_3UTR_length)

all_data_tsg$hg38_5UTR_length_log=log10(genes_infor[all_data_tsg$ensembl_gene_id,]$hg38_5UTR_length)
all_data_tsg$cds_length_log=genes_infor[all_data_tsg$ensembl_gene_id,]$cds_length_log

all_data_tsg$rbp_site_binding_density2=genes_infor[all_data_tsg$ensembl_gene_id,]$rbp_site_binding_density2
all_data_tsg$rbp_site_binding_density2_support=genes_infor[all_data_tsg$ensembl_gene_id,]$rbp_site_binding_density2_support

exp_norm_TSMI_TCGA<-aggregate(all_data_tsg$norm_exp,list(all_data_tsg$ensembl_gene_id),get_TSEI_log2)
rownames(exp_norm_TSMI_TCGA)<-exp_norm_TSMI_TCGA[,1]
genes_infor$TSMI_tissue_specificity_TCGA=exp_norm_TSMI_TCGA[genes_infor$ensembl_gene_id,2]

all_data_tsg_all_genes=all_data_tsg
all_data_tsg<-all_data_tsg[which(all_data_tsg$Gene.Type!="Oncogene"),]


all_data_tsg$motif_count_density<-genes_infor[all_data_tsg$ensembl_gene_id,"motif_count_density"]
all_data_tsg_all_genes$motif_count_density<-genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"motif_count_density"]

all_data_tsg_all_genes$motif_count<-log10(genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"motif_count"])
genes_infor$rbp_motif<-paste(genes_infor$rbp,genes_infor$motif)





filter_low_data<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1
    #print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="TSG" & temp_data$diff_exp<0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}

all_data_tsg_filtered=filter_low_data(all_data_tsg )
all_data_tsg_filtered$diff_exp=abs(all_data_tsg_filtered$diff_exp)


save.image("./data/combined_gene_infor_expression_data.RData")

