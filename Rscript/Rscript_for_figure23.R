load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')
library(dplyr)
rownames(gene_age_infor)<-gene_age_infor$Age.class
genes_infor$age_type3<-as.character(genes_infor$age_type2)
genes_infor$age_type3[genes_infor$age_type3%in% c("(67.6 - 355.7]","(0 - 67.6]")]="(0 - 355.7]"
p1_3UTR_age<-plotAgelength(genes_infor,"hg38_3UTR_length_log","Log10(Length of 3'UTR)",0,"age_type3")
genes_infor$rbp<-ifelse(is.na(genes_infor$rbp_site_binding_density2),"RBP-","RBP+")
genes_infor$miRNA<-ifelse(is.na(genes_infor$mirna_binding_3utr_density),"miRNA-","miRNA+")

genes_infor=genes_infor[genes_infor$Gene.Type!="Oncogene",]

p.ff_tcga_length<-plotForest(all_data_tsg_all_genes,"norm_exp","hg38_3UTR_length_log","Pearson's r(Log10 3'UTR Length VS\nGene Expression) in Normal samples(TCGA)",0.5,1.1,"1.5cm")

age_levels<-c("<67.6",gene_age_infor[8:26,4])
all_data_tsg$age_time<-as.numeric(gene_age_infor[as.character(all_data_tsg$Age_class),4])
all_data_tsg$age_time[all_data_tsg$age_time<=67.6]="<67.6"
age_cor_len<-get_cor_for_each_age(all_data_tsg,"hg38_3UTR_length_log","norm_exp","age_time","pearson")
age_cor_len$age<-factor(age_cor_len$age,levels = age_levels)

p_cor_age_len<-ggplot(data=age_cor_len ,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25)+theme_hd()+xlab("Gene age class(million years)")+theme(axis.text.x =element_text(  angle=90,hjust=1))+ylab("Pearson's r(Log10(3'UTR length) VS\n Gene Expression)")+ggtitle("Normal samples\nTCGA RNA-seq datasets")




p_berbp_cor<-plot_scatter_6(genes_infor[genes_infor$Gene.Type!="Oncogene",],"hg38_3UTR_length_log","rbp_site_binding2","Log10(3'UTR length)","Count of RBPs \nbinding sits in 3'UTR",0.3,1,"black")+theme(plot.title = element_text(hjust = 0.5)) +theme_hd()




p_berbp_density=ggplot(data = genes_infor[genes_infor$Gene.Type!="Oncogene",],aes(y=rbp_site_binding_density2,x=Gene.Type,fill=Gene.Type))+geom_boxplot(width=0.4)+ylab("Log10(RBPs binding \ndensity in 3'UTR)")+geom_signif(comparisons = list(c("TSG","Non-Cancer")),map_signif_level=T,col="black", test='wilcox.test',fontface="bold",y_position = -0.4)+theme_hd()+theme(axis.title.x = element_blank(),legend.position = "none" )+scale_fill_manual(values=c("#1874CD","#EE2C2C"))
#axis.text.x=element_text(angle = 30,hjust = 1),

p_berbp_count=ggplot(data = genes_infor[genes_infor$Gene.Type!="Oncogene",],aes(y=rbp_site_binding2,x=Gene.Type,fill=Gene.Type))+geom_boxplot(aes(col=I("gray90")),width=0.4)+ylab("No.of RBPs binding site in 3'UTR)")+geom_signif(comparisons = list(c("TSG","Non-Cancer")),map_signif_level=T,col="black", test='wilcox.test',fontface="bold")+ xlab("Dataset")+theme_hd()+theme(axis.text.x=element_text(colour = "white") )+xlab("")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+xlab("Gene Types")

p_berbp_typecount=ggplot(data = genes_infor[genes_infor$Gene.Type!="Oncogene",],aes(y=rbp_type_count,x=Gene.Type,fill=Gene.Type))+geom_boxplot(aes(col=I("gray90")),width=0.4)+ylab("No.of RBPs  in 3'UTR)")+geom_signif(comparisons = list(c("TSG","Non-Cancer")),map_signif_level=T,col="black", test='wilcox.test',fontface="bold")+ xlab("Dataset")+theme_hd()+theme(axis.text.x=element_text(colour = "white") )+xlab("")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+xlab("Gene Types")


p_berbp_3utrlen=ggplot(data = genes_infor[genes_infor$Gene.Type!="Oncogene",],aes(y=hg38_3UTR_length_log,x=rbp,fill=rbp))+geom_boxplot(width=0.4)+ylab("log10 3'UTR length")+geom_signif(comparisons = list(c("RBP+","RBP-")),map_signif_level=T,col="black", test='wilcox.test',fontface="bold")+theme_hd()+xlab("RBPs in 3'UTR")+scale_fill_manual("", values=c("#76EE00","darkgreen"))+theme(legend.position = "none")






all_data_tsg_all_genes$rbp_site_binding2=log10(genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"rbp_site_binding2"])
all_data_tsg_all_genes$rbp_type_count=genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"rbp_type_count"]
p.ff_tcga_berbp<-plotForest(all_data_tsg_all_genes,"norm_exp","rbp_site_binding_density2","Pearson's r(Log10(RBPs binding density in 3'UTR)\n VS Gene Expression) in Normal samples(TCGA )",0.5,1.1,"1cm")
p.ff_tcga_berbp2<-plotForest(all_data_tsg_all_genes,"norm_exp","rbp_site_binding2","Pearson's r(Log10(No of RBPs binding site in 3'UTR)\n VS Gene Expression) in Normal samples(TCGA )",0.5,0.9,"1cm")
p.ff_tcga_berbp3<-plotForest(all_data_tsg_all_genes,"norm_exp","rbp_type_count","Pearson's r(No.of RBPs  in 3'UTR)\n VS Gene Expression) in Normal samples(TCGA )",0.5,0.9,"1cm")



img <- readPNG("./data/meme_RNA_3/logo1.png")

g <- rasterGrob(img, interpolate=TRUE)
s0=qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=0, xmax=10, ymin=5, ymax=10) +
  geom_point(aes(col=I("white")))+theme_void()


s0=s0+annotate("text",x=2.4,y=4,hjust=.2,size=4,label="Frequency in          Frequency in    P-value    ",fontface="bold")
s0=s0+annotate("text",x=1.7,y=3,hjust=.2,size=4,label="Non-cancer            TSG",fontface="bold")
s0=s0+annotate("text",x=2.5,y=2,hjust=.2,size=4,label="42%                          55%                  = 8.64e-15",fontface="bold")


fimo <- read.delim("./data/meme_RNA_3/fimo.tsv", comment.char="#")
geneid=substr(fimo$sequence_name,17,31)
genes_infor$motif=rep("C-rich_motif-",length(genes_infor$ensembl_gene_id))
genes_infor$motif[genes_infor$ensembl_gene_id %in% geneid]="C-rich_motif+"
temp_x<-table(geneid)
genes_infor$motif_count=as.numeric(temp_x[genes_infor$ensembl_gene_id])
genes_infor$motif_count_density=log10(genes_infor$motif_count/genes_infor$hg38_3UTR_length)
genes_infor1=genes_infor[which(genes_infor$Gene.Type !="Oncogene"),]
all_data_tsg$motif_count_density<-genes_infor[all_data_tsg$ensembl_gene_id,"motif_count_density"]
all_data_tsg_all_genes$motif_count_density<-genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"motif_count_density"]

all_data_tsg_all_genes$motif_count<-log10(genes_infor[all_data_tsg_all_genes$ensembl_gene_id,"motif_count"])
genes_infor$rbp_motif<-paste(genes_infor$rbp,genes_infor$motif)





x1=table(genes_infor1$motif , genes_infor1$Gene.Type)
fisher.test(x1)
x1[2,]/apply(x1,2,sum)


p_ff_motif<-plotForest(all_data_tsg_all_genes,"norm_exp","motif_count_density","Pearson's r(Log10(C-rich motif density in 3'UTR)\n VS Gene Expression) in Normal samples(TCGA )",0.5,1.1,"1cm")


plot_scatter_2(all_data_tsg_all_genes[!is.na(all_data_tsg_all_genes$age_type2) & all_data_tsg_all_genes$cancer=="COAD",],"motif_count_density","norm_exp","Log10(C-rich motif density in 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA-COAD normal ")+theme_hd()+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~age_type2,nrow=1)

p_count_cor<-plot_scatter_6(genes_infor,"hg38_3UTR_length_log","motif_count","Log10(Length 0f 3'UTR)","C-rich motif\n count per gene",0.3,1,"black")+theme(plot.title = element_text(hjust = 0.5)) +theme_hd()

all_data_tsg$motif<-genes_infor[as.character(all_data_tsg$ensembl_gene_id),]$motif
all_data_tsg$rbp<-genes_infor[as.character(all_data_tsg$ensembl_gene_id),]$rbp
all_data_tsg$rbp_motif<-genes_infor[as.character(all_data_tsg$ensembl_gene_id),]$rbp_motif
countFunction <- function(x){
  return(data.frame(y=rep(0,length(x)),label=round(length(x),2)))}

comp=list(c("C-rich_motif-","C-rich_motif+"))
ggplot(data = all_data_tsg,aes(y=norm_exp,x=motif,fill=motif))+geom_boxplot()+geom_signif(comparisons = comp,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Log2(TPM+1) expression ")+facet_wrap(~cancer,nrow = 1)
sig=c(rep("***",13))


comp=list(c("C-rich_motif-","C-rich_motif+"))
s1_motif=ggplot(data = all_data_tsg,aes(y=norm_exp,x=cancer,fill=motif))+geom_boxplot(aes(col=I("gray90")))+scale_fill_manual("Motif in 3'UTR", values=c("#76EE00","darkgreen"))+ coord_flip()+annotate("text",1:13,15,label=sig,col="red")+theme(legend.position = "bottom")+ylab("Gene Expression")+ggtitle("Normal samples")+ xlab("Dataset")+theme_hd()


comp=list(c("RBP-","RBP+"))
all_data_tsg$rbp=ifelse(is.na(all_data_tsg$rbp_site_binding_density2),"RBP-","RBP+")
ggplot(data = all_data_tsg,aes(y=norm_exp,x=rbp,fill=rbp))+geom_boxplot()+geom_signif(comparisons = comp,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1, test.args = list(alternative = "less"))+ylab("Log2(TPM+1) expression ")+facet_wrap(~cancer,nrow = 1)
s1_rbp=ggplot(data = all_data_tsg,aes(y=norm_exp,x=cancer,fill=rbp))+geom_boxplot(aes(col=I("gray90")))+scale_fill_manual("RBPs in 3'UTR", values=c("#76EE00","darkgreen"))+ coord_flip()+annotate("text",1:13,15,label=sig,col="red")+theme(legend.position = "bottom")+ylab("Gene Expression")+ggtitle("Normal samples")+ xlab("Dataset")+theme_hd()


genes_infor=genes_infor[genes_infor$Gene.Type!="Oncogene",]
x1=table(genes_infor$rbp , genes_infor$Gene.Type)
fisher.test(x1)
x1[2,]/apply(x1,2,sum)

x1=table(genes_infor$miRNA , genes_infor$Gene.Type)
fisher.test(x1)
x1[1,]/apply(x1,2,sum)


p0<-ggarrange(p.ff_tcga_length,p_cor_age_len, p_berbp_3utrlen,nrow=1,ncol=3,labels=c(letters[c(1:3)]), font.label = list(size = 18, color = "black", face = "bold"),widths = c(2,1.3,0.7))
p01<-ggarrange(p_berbp_cor,p_berbp_density,labels=c(letters[c(5:6)]),nrow=2,ncol=1, font.label = list(size = 18, color = "black", face = "bold"))
p02<-ggarrange(s1_rbp,p01,nrow=1,ncol=2,labels=c(letters[4],""),font.label = list(size = 18, color = "black", face = "bold"))
p1<-ggarrange(p02,p.ff_tcga_berbp,nrow=1,ncol=2,labels=c("",letters[7]),font.label = list(size = 18, color = "black", face = "bold"))

p12<-ggarrange(s0,p_count_cor, nrow=2,ncol=1,labels=c(letters[8:9]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1),common.legend = T,legend="bottom")
p2<-ggarrange(p12,s1_motif,p_ff_motif,nrow=1,ncol=3,labels=c("",letters[10],letters[11]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1,2))

pall<-ggarrange(p0,p1,p2,nrow=3,ncol=1)


tiff("./result/figure_2_3UTR.tiff",width = 16.5,height = 16,res=300,units="in",compression = "lzw")
print(pall)
dev.off()



###################################################################figure3

rbp_count=as.data.frame(table(genes_infor[which(genes_infor$rbp=="RBP+"),]$age_type2,genes_infor[which(genes_infor$rbp=="RBP+"),]$Gene.Type))
rbp_all=as.data.frame(table(genes_infor$age_type2,genes_infor$Gene.Type))
merge_rbp=merge(rbp_count,rbp_all,by=c("Var1","Var2"))
merge_rbp$ratio=merge_rbp$Freq.x/merge_rbp$Freq.y

p_rbp<-ggplot(data=merge_rbp[merge_rbp$Var2!="Oncogene",],aes(x=Var1,y=ratio,fill=Var2))+geom_bar(stat = "identity", position=position_dodge())+scale_fill_manual("Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab("Gene Age Groups(million years)")+ylab("Ratio of genes with RBP binding site ")+theme_hd()+theme(axis.text.x = element_text(angle = 90))
comp=list(c("RBP-","RBP+"))

all_data_tsg$rbp<-factor(all_data_tsg$rbp,levels = c("RBP-","RBP+"))
p_exp_rbp<-plotSubSig(all_data_tsg[which(all_data_tsg$Gene.Type!="Oncogene" & !is.na(all_data_tsg$age_type2) & all_data_tsg$cancer=="COAD"),],"age_type2", "norm_exp","rbp","", cols1 =c("#76EE00","darkgreen"),"Gene Age Groups(million years)","Gene Expression") + labs(fill = "RBPs in 3'UTR")+ggtitle("TCGA-COAD normal")+theme(axis.text.x = element_text(angle = 20))


p_age_motif<-ggplot(genes_infor[!is.na(genes_infor$age_type3),],aes(x=age_type3,y=motif_count_density,fill=age_type3))+geom_boxplot(show.legend = FALSE)+geom_signif(comparisons=getComp(unique(genes_infor$age_type3[!is.na(genes_infor$age_type2)])), test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +theme_hd()+scale_fill_manual(values=brewer.pal(n = 5, name = 'Blues')[2:5])+theme(axis.text.x = element_text(angle = 90,hjust = 1))+ylab("Log10(C-rich motif density in 3'UTR)")+xlab("Gene Age Groups\n(million years)")

p_exp<-plotSubSig(all_data_tsg[which(all_data_tsg$Gene.Type!="Oncogene" & !is.na(all_data_tsg$age_type2) & all_data_tsg$cancer=="COAD"),],"age_type2", "norm_exp","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age Groups(million years)","Gene Expression") + labs(fill = "Gene.Type")+ggtitle("TCGA-COAD normal ")+theme(legend.position = "bottom")

pcor1<-plot_scatter_2(all_data_tsg_all_genes[!is.na(all_data_tsg_all_genes$age_type2) & all_data_tsg_all_genes$cancer=="COAD",],"hg38_3UTR_length_log","norm_exp","Log10(Length of 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA-COAD normal ")+theme_hd()+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~age_type2,nrow=1)

p_dif_rbp<-plotSubSigOneSide(genes_infor[!is.na(genes_infor$age_type2) & genes_infor$Gene.Type!="Oncogene",],"age_type2", "rbp_site_binding_density2","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age Groups (million years)","Log10(RBPs binding density in 3'UTR)",side="less") + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")

p_dif_rbp_count<-plotSubSigOneSide(genes_infor[!is.na(genes_infor$age_type2) & genes_infor$Gene.Type!="Oncogene",],"age_type2", "rbp_site_binding2","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age Groups (million years)","Log10(RBPs binding site in 3'UTR)",side="less") + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")


pcor_rerbp<-plot_scatter_2(all_data_tsg_all_genes[!is.na(all_data_tsg_all_genes$age_type2) & all_data_tsg_all_genes$cancer=="COAD",],"rbp_site_binding_density2","norm_exp","Log10(RBPs binding density in 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA-COAD normal ")+theme_hd()+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~age_type2,nrow=1)


motif_count=as.data.frame(table(genes_infor[which(genes_infor$motif=="C-rich_motif+"),]$age_type2,genes_infor[which(genes_infor$motif=="C-rich_motif+"),]$Gene.Type))
motif_all=as.data.frame(table(genes_infor$age_type2,genes_infor$Gene.Type))
merge_motif=merge(motif_count,motif_all,by=c("Var1","Var2"))

merge_motif$ratio=merge_motif$Freq.x/merge_motif$Freq.y
p_motif<-ggplot(data=merge_motif[merge_motif$Var2!="Oncogene",],aes(x=Var1,y=ratio,fill=Var2))+geom_bar(stat = "identity", position=position_dodge())+scale_fill_manual("Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab("Gene Age Groups(million years)")+ylab("Ratio of genes with C-rich Motif ")+theme_hd()+theme(axis.text.x = element_text(angle = 90))


motif_ratio=as.data.frame(genes_infor[genes_infor$Gene.Type!="Oncogene" & !is.na(genes_infor$age_type2),] %>% dplyr::count(motif,age_type2,Gene.Type) %>% dplyr::group_by(age_type2,Gene.Type) %>%dplyr::mutate(freq = n / sum(n)))

p_motif<-ggplot(data=motif_ratio[motif_ratio$motif=="C-rich_motif+",],aes(x=age_type2,y=freq,fill=Gene.Type))+geom_bar(stat = "identity", position=position_dodge())+scale_fill_manual("Gene.Type",values=c("#1874CD","#EE2C2C"))+xlab("Gene Age Groups(million years)")+ylab("Ratio of genes with C-rich Motif ")+theme_hd()+theme(axis.text.x = element_text(angle = 90))







comp=list(c("C-rich_motif-","C-rich_motif+"))
p_exp_motif<-plotSubSig(all_data_tsg[which(all_data_tsg$Gene.Type!="Oncogene" & !is.na(all_data_tsg$age_type2) & all_data_tsg$cancer=="COAD"),],"age_type2", "norm_exp","motif","", cols1 =c("#76EE00","darkgreen"),"Gene Age Groups(million years)","Gene Expression") + labs(fill = "Motif in 3'UTR")+ggtitle("TCGA-COAD normal")+theme(axis.text.x = element_text(angle = 20))

pcor_motif<-plot_scatter_2(all_data_tsg_all_genes[!is.na(all_data_tsg_all_genes$age_type2) & all_data_tsg_all_genes$cancer=="COAD",],"motif_count_density","norm_exp","Log10(C-rich motif density in 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA-COAD normal ")+theme_hd()+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~age_type2,nrow=1)




tcgs_cor_age_tsg_rbp<-get_cor_agetype(all_data_tsg,"norm_exp","rbp_site_binding_density2","age_type2")

p_cor_age_tcgs_rbp<-ggplot(data=tcgs_cor_age_tsg_rbp,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10(RBPs binding density\nin 3'UTR) VS Gene Expression)")+ggtitle("Normal samples\nTCGA RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))

cor_age_tsg<-get_cor_agetype(all_data_tsg_all_genes,"norm_exp","motif_count_density","age_type2")

p_cor_age_tsg_motif<-ggplot(data=cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10(C-rich motif density\nin 3'UTR) VS Gene Expression)")+ggtitle("Normal samples\nTCGA RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))



cor_age_tsg<-get_cor_agetype(all_data_tsg_all_genes,"norm_exp","rbp_type_count","age_type2","kendall")

p_cor_age_tsg_rbp_type<-ggplot(data=cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Kendall's tau(No.of RBPs \nin 3'UTR) VS Gene Expression)")+ggtitle("Normal samples\nTCGA RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))




pm<-ggarrange(p_motif,p_cor_age_tsg_motif,nrow=1,ncol=2,labels=c(letters[7:8]), font.label = list(size = 18, color = "black", face = "bold"))
p0<-ggarrange(p_rbp,p_cor_age_tcgs_rbp,labels=c(letters[4:5]),nrow=1,ncol=2, font.label = list(size = 18, color = "black", face = "bold"))

p0<-ggarrange(p_exp+theme(axis.text.x = element_text(angle = 20)),pcor1,labels=c(letters[1:2]),nrow=1,ncol=2, font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,2.2))

p2=ggarrange(p_exp_rbp,p_rbp,p_exp_motif,p_motif,nrow = 2,ncol=2,labels=c(letters[3:6]), font.label = list(size = 18, color = "black", face = "bold"))

p3<-ggarrange(p0,p2,nrow=2,heights = c(1,2))
tiff("./result/figure_3_3UTR.tiff",width = 15,height = 17,res=300,units="in",compression = "lzw")
print(p3)
dev.off()








#supplymentary figures

p_berbp_3utr= plotSubSig(genes_infor,"rbp_motif", "hg38_3UTR_length_log","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene groups","Log10(Length of 3'UTR)") + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 78,hjust = 1))

p_berbp_3utr2= ggplot(data = genes_infor,aes(y=hg38_3UTR_length_log,x=rbp_motif,fill=rbp_motif))+geom_boxplot(width=0.4)+ylab("Log10(Length of 3'UTR)")+geom_signif(comparisons = getComp(unique(all_data_tsg$rbp_motif)),map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,fontface="bold")+theme_hd()+xlab("RBPs in 3'UTR")+theme_hd()+xlab("Gene Groups")+theme(axis.text.x = element_blank(),legend.position = "none")


rbp_motif_ratio=as.data.frame(genes_infor[genes_infor$Gene.Type!="Oncogene",] %>% count(Gene.Type,rbp_motif) %>% group_by(Gene.Type) %>%mutate(freq = n / sum(n)))

p_rbp_motif<-ggplot(data=rbp_motif_ratio,aes(x=Gene.Type,y=freq,fill=rbp_motif))+geom_bar(stat = "identity", position=position_dodge())+theme_hd()+ylab("Ratio of genes ")+theme(legend.position = "none")




##############################################################################
tiff("./result/supplymentary_figure3.tiff",width = 10,height = 10,res=300,units="in",compression = "lzw")
sp1=plot_scatter_7(all_data_tsg_all_genes[all_data_tsg_all_genes$cancer!="COAD",],"hg38_3UTR_length_log","norm_exp","Log10(Length of 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA RNA-seq data Normal samples")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4)
print(sp1)
dev.off()


tiff("./result/supplymentary_figure4.tiff",width = 10,height = 10,res=300,units="in",compression = "lzw")
sp1=plot_scatter_7(all_data_tsg_all_genes[all_data_tsg_all_genes$cancer!="COAD",],"rbp_site_binding_density2","norm_exp","Log10(RBPs binding density in 3'UTR)", "Gene Expression",0.01,15,"black")+ggtitle("TCGA RNA-seq data Normal samples")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4)
print(sp1)
dev.off()




g1 <- rasterGrob(readPNG("/media/huangdan/hardisk0/HD/my_data/HD_evolution/result_UTR/meme_RNA_3/logo_PCBP2.png"), interpolate=TRUE)
g2 <- rasterGrob(readPNG("/media/huangdan/hardisk0/HD/my_data/HD_evolution/result_UTR/meme_RNA_3/logo_HRB87F.png"), interpolate=TRUE)
g3 <- rasterGrob(readPNG("/media/huangdan/hardisk0/HD/my_data/HD_evolution/result_UTR/meme_RNA_3/logo_HNRNPH2.png"), interpolate=TRUE)


s01=qplot(1:10, 1:10, geom="blank") +annotation_custom(g1, xmin=-1, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(aes(col=I("white")))+theme_void()+annotate("text",x=7,y=10,hjust=.2,size=4,label="(PCBP2) q-value:  2.01e-02")
s02=qplot(1:10, 1:10, geom="blank") +annotation_custom(g2, xmin=-1, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(aes(col=I("white")))+theme_void()+annotate("text",x=7,y=10,hjust=.2,size=4,label="(HRB87F) q-value:3.54e-02")
s03=qplot(1:10, 1:10, geom="blank") +annotation_custom(g3, xmin=-1, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(aes(col=I("white")))+theme_void()+annotate("text",x=7,y=10,hjust=.2,size=4,label="(HNRNPH2) q-value:  3.59e-04")


p_all0_supply<-ggarrange(s03,s01,s02 ,nrow  = 3,ncol=1)
tiff("./result/supplymentary_figure5.tiff",width = 5,height = 7,res=300,units="in",compression = "lzw")
print(p_all0_supply)
dev.off()



tiff("./result/supplymentary_figure7.tiff",width = 13,height = 16,res=300,units="in",compression = "lzw")
p<-plotSubSigall4(all_data_tsg[all_data_tsg$cancer!="COAD" & !is.na(all_data_tsg$age_type2),],"age_type2", "norm_exp","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age Groups(million years)","Gene Expression") 
print(p)
dev.off()



tcgs_cor_age_tsg<-get_cor_agetype(all_data_tsg,"norm_exp","hg38_3UTR_length_log","age_type2")

p_cor_age_tcgs<-ggplot(data=tcgs_cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age Groups(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10 3'UTR Length VS\n Gene Expression)")+ggtitle("Normal samples\nTCGA RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))

all_data_tsg$age_type3<-as.character(all_data_tsg$age_type2)
all_data_tsg$age_type3[all_data_tsg$age_type3%in% c("(67.6 - 355.7]","(0 - 67.6]")]="(0 - 355.7]"
all_data_tsg$age_type3<-as.factor(all_data_tsg$age_type3)
#,pcor_rerbp,pcor_motif,

tiff("./result/supplymentary_figure8.tiff",width = 7,height = 7,res=300,units="in",compression = "lzw")
print(p_cor_age_tcgs)
dev.off()

tiff("./result/supplymentary_figure9.tiff",width = 13,height = 16,res=300,units="in",compression = "lzw")
p<-plotSubSigall4(all_data_tsg[all_data_tsg$cancer!="COAD" & !is.na(all_data_tsg$age_type2),],"age_type2", "norm_exp","rbp","", cols1 = c("#76EE00","darkgreen"),"Gene Age Groups(million years)","Gene Expression") 
print(p)
dev.off()

tiff("./result/supplymentary_figure10.tiff",width = 13,height = 16,res=300,units="in",compression = "lzw")
p<-plotSubSigall4(all_data_tsg[all_data_tsg$cancer!="COAD" & !is.na(all_data_tsg$age_type2),],"age_type2", "norm_exp","motif","", cols1 = c("#76EE00","darkgreen"),"Gene Age Groups(million years)","Gene Expression") 
print(p)
dev.off()


