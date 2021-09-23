
load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')
filter_low_data_Oncogene<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1
    #print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="Oncogene" & temp_data$diff_exp>0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}


all_data_tsg_filtered_Onco=filter_low_data_Oncogene(all_data_tsg_all_genes)

p.ff_diff_onco<-plotForest(all_data_tsg_filtered_Onco,"diff_exp","hg38_3UTR_length_log","Pearson's r(Log10 3'UTR Length VS log2FC)\nof Up-regulated Oncogenes in TCGA",0.5,1.5,"1cm")





filter_low_dataOncoArray<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=(quantile(temp_data$norm_exp,0.2,na.rm=T)+quantile(temp_data$tumor_exp,0.2,na.rm=T))/2
    print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="Oncogene" & temp_data$diff_exp>0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}

all_data_other_cohort_filteredOnCo=filter_low_dataOncoArray(all_data_other_cohort  )



p.ff_diffveriOnco<-plotForest(all_data_other_cohort_filteredOnCo,"diff_exp","hg38_3UTR_length_log","Pearson's r(Log10 3'UTR Length VS log2FC)\nof Up-regulated Oncogenes in validation datasets",0.9,1.5,"1cm")

p0_3utr_all<-plot_all_genes(genes_infor,"hg38_3UTR_length_log","Log10(Length of 3'UTR)") 
p0_3utr_exp<-plot_all_genes(all_data_tsg_all_genes[which( all_data_tsg_all_genes$cancer=="LUSC"),],"norm_exp","Gene Expression") 

p_num<-ggplot(genes_infor[genes_infor$Gene.Type=="Oncogene" & !is.na(genes_infor$age_type2),],aes(x=age_type2))+geom_bar(position = "dodge")+theme_hd()+theme(axis.text.x = element_text(angle = 30,hjust = 1),strip.background = element_blank())+ylab("No. of Oncogenes")+xlab("Gene Age Groups(million years)")










tiff("./result/supplymentary_figure16.tiff",width = 17,height = 14,res=300,units="in",compression = "lzw")
p0=ggarrange( p0_3utr_all,p.ff_diffveriOnco,nrow=1,ncol=2,labels = c(letters[1:2]), font.label = list(size = 18, color = "black", face = "bold"))
p1=ggarrange( p.ff_diff_onco, p_num,nrow=1,ncol=2,labels = c(letters[3:4]), font.label = list(size = 18, color = "black", face = "bold"))

ggarrange(p0,p1,heights = c(2,2),nrow=2,ncol=1)
dev.off()



tiff("./result/supplymentary_figure15.tiff",width = 13,height = 13,res=300,units="in",compression = "lzw")
sp1=plot_scatter_2(all_data_tsg_all_genes,"mirna_binding_3utr_density","norm_exp","Log10(microRNA binding density in 3'UTR)",
                   "Gene Expression",0.01,15,"black")+ggtitle("TCGA RNA-seq data Normal samples")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4)
print(sp1)

dev.off()





