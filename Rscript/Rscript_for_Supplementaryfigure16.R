
load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')


SigMut_List <- read.csv("./data/SigMut_List.csv")
SigMut_List$sig<--log10(SigMut_List$p.val)
SigMut_List$Gene.Type<-genes_infor[SigMut_List$Ensembl.ID,]$Gene.Type
SigMut_List$hg38_3UTR_length_log<-genes_infor[SigMut_List$Ensembl.ID,]$hg38_3UTR_length_log


p.ff_tcga_sig3<-plotForest(SigMut_List_TSG_all[SigMut_List_TSG_all$Gene.Type.x=="TSG",],"sig","hg38_3UTR_length_log.y","Pearson's r(Log10 3'UTR Length VS\n-Log10(MutsigCV_Pvalue)) in Normal samples(TCGA)",0.5,1.1,"1.5cm")



sm2=plot_scatter_2(genes_infor,"hg38_3UTR_length_log","cds_length_log","Log10(3'UTR length)","Log10(CDS length)",0.4,1,"black")+theme(plot.title = element_text(hjust = 0.5,face="bold"))+ theme(strip.text = element_text(size=5),axis.title =element_text(size=14))
print(sm2)

sm3=plot_scatter_2(genes_infor[genes_infor$Gene.Type=="TSG",],"hg38_3UTR_length_log","cds_length_log","Log10(3'UTR length)","Log10(CDS length)",0.4,1,"black")+theme(plot.title = element_text(hjust = 0.5,face="bold"))+ theme(strip.text = element_text(size=5),axis.title =element_text(size=14))+ggtitle("TSG")
print(sm3)



x=table(all_data_tsg_filtered$ensembl_gene_id)
genes_tsg<-genes_infor[which(genes_infor$ensembl_gene_id %in% names(x) & !is.na(genes_infor$hg38_3UTR_length_log)),]
genes_tsg$cancer_speci<-x[genes_tsg$ensembl_gene_id]
genes_tsg$cancer_speci_ratio<-genes_tsg$cancer_speci/13




sd1=plot_scatter_2(genes_tsg,"hg38_3UTR_length_log","cancer_speci_ratio","Log10(Length of 3'UTR)", "Cancer Heterogeneity",0,4.5,"black")+ggtitle(paste("Down-regulated TSGs"))+theme(plot.title = element_text(hjust = 0.5)) +theme_hd()+ylim(0,1.1)
print(sd1)


p.ff_diff<-plotForest(all_data_tsg_filtered[all_data_tsg_filtered$ensembl_gene_id %in% names(x[x>6]),],"diff_exp","hg38_3UTR_length_log","Pearson's r(Log10_3UTR_Length VS abs(log2FC))\n Universally Down-regulated TSG in TCGA",0.9,1,"5cm")



library(readxl )
percentage_cancer <-as.data.frame( read_excel("./data/percentage_Comprehensive Characterization of Cancer Driver Genes and Mutations.xlsx"))
percentage_cancer$cancer[percentage_cancer$cancer=="COADREAD"]="COAD"
percentage_cancer=percentage_cancer[-which(percentage_cancer$cancer=="STAD"),]
rownames(percentage_cancer)<-percentage_cancer$cancer
percentage_cancer$type="Unknown"
percentage_cancer$type[as.numeric(percentage_cancer$TSG)>=0.60]="TSG driver cancers"
percentage_cancer$type[as.numeric(percentage_cancer$Oncogene)>=0.60]="Oncogene driver cancers"
percentage_cancer$cancer=factor(percentage_cancer$cancer,levels = percentage_cancer$cancer[order(percentage_cancer$TSG)])


cor_plot_len_diff<-get_cor_for_each_cancer_all(all_data_tsg_filtered[all_data_tsg_filtered$ensembl_gene_id %in% names(x[x>6]),],"diff_exp","hg38_3UTR_length_log","pearson")
all_cor=cbind(type="3utr_diff",cor_plot_len_diff)

all_cor$percent_cancer_driver<-percentage_cancer[as.character( all_cor$Cancer),2]

all_cor$cancer_type<-percentage_cancer[as.character( all_cor$Cancer),"type"]
all_cor$r2<-all_cor$cor^2


pd1 <- ggscatter(all_cor[all_cor$type=="3utr_diff",] ,  y = "cor", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's r(Log10 3'UTR Length\nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Universally Down-regulated TSG")+ylim(-0.35,0)
#+ylim(0,0.1)+xlim(0.15,0.9)



p1=ggarrange(sm2,p.ff_tcga_sig3,sd1,nrow = 1,ncol=3,labels = letters[c(1:3)],widths = c(0.5,1.5,1),font.label = list(size = 18, color = "black", face = "bold"))
p2=ggarrange(p.ff_diff,pd1, nrow = 1,ncol=2,labels = letters[c(4,5)],widths = c(2,1),font.label = list(size = 18, color = "black", face = "bold"))

tiff("./result/supplymentary_figure16.tiff",width = 24,height = 10,res=300,units="in",compression = "lzw")
print(ggarrange(p1,p2,nrow=2,ncol=1))
dev.off()










