load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')

filter_low_data<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=(quantile(temp_data$norm_exp,0.2,na.rm=T)+quantile(temp_data$tumor_exp,0.2,na.rm=T))/2
    print(quantile_cutoff)
    temp_data=temp_data[which(temp_data$Gene.Type=="TSG" & temp_data$diff_exp<0  & !(temp_data$norm_exp < quantile_cutoff & temp_data$tumor_exp < quantile_cutoff )),]  
    result=rbind(result,temp_data)
  }
  
  return(result)
}


all_data_other_cohort_filtered=filter_low_data(all_data_other_cohort  )
all_data_other_cohort_filtered$diff_exp<-abs(all_data_other_cohort_filtered$diff_exp)




all_data_other_cohort_filtered$rbp_site_binding2=log10(genes_infor[all_data_other_cohort_filtered$ensembl_gene_id,"rbp_site_binding2"])
all_data_other_cohort_filtered$rbp_type_count=log10(genes_infor[all_data_other_cohort_filtered$ensembl_gene_id,"rbp_type_count"])



pcor1=plotForest(all_data_other_cohort_filtered,"diff_exp","rbp_site_binding2","Pearson's r(Log10 No.of RBP binding site VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")
pcor2=plotForest(all_data_other_cohort_filtered,"diff_exp","rbp_type_count","Pearson's r(Log10 No.of RBPs VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")
pcor3=plotForest(all_data_other_cohort_filtered,"diff_exp","rbp_site_binding_density2","Pearson's r(Log10 rbp_site_binding_density VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")


genes_infor$rbp<-ifelse(is.na(genes_infor$rbp_site_binding_density2),"RBP-","RBP+")
all_data_other_cohort_filtered$rbp=genes_infor[all_data_other_cohort_filtered$ensembl_gene_id,"rbp"]
all_data_other_cohort_filtered$motif=genes_infor[all_data_other_cohort_filtered$ensembl_gene_id,"motif"]

p1=ggplot(data=all_data_other_cohort_filtered,aes(x=rbp,y=diff_exp,fill=rbp))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("With","Without")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold")+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+ theme(plot.title = element_text(face = "bold"),axis.text.x = element_text(angle = 90))+xlab("RBP binding site in 3'UTR")+facet_wrap(~cancer,nrow=1)

sig_rbp<-get_diff(all_data_other_cohort_filtered,"rbp","diff_exp","RBP+")
pff_rbp_all<-plotForesttwogroup(sig_rbp,"Abs(Log2FC) of Down-regulated TSG in TCGA","1mm","RBP+","RBP-")
pff_rbp_in<-pff_rbp_all$p1

sig=c(rep("***",53))
comp=list(c("RBP-","RBP+"))
GTEXmedian_tpm$rbp=ifelse(is.na(GTEXmedian_tpm$rbp_site_binding_density2),"RBP-","RBP+")
ggplot(data = GTEXmedian_tpm,aes(y=norm_exp,x=rbp,fill=rbp))+geom_boxplot(aes(col=I("gray")))+geom_signif(comparisons = comp,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,y_position = 10)+ylab("Log2(TPM+1) expression ")+facet_wrap(~cancer)
s1_rbp=ggplot(data = GTEXmedian_tpm,aes(y=norm_exp,x=cancer,fill=rbp))+geom_boxplot(aes(col=I("gray90")))+scale_fill_manual("RBPs in 3'UTR", values=c("#76EE00","darkgreen"))+ coord_flip()+annotate("text",1:53,15,label=sig,col="red")+theme(legend.position = "bottom")+ylab("Gene Expression")+ggtitle("Normal samples")+ xlab("Dataset")+theme_hd()







tiff("./result/test_supfigure_4_3UTR.tiff",width = 20,height = 12,res=300,units="in",compression = "lzw")
print(ggarrange(pcor1,pcor2,pcor3,p1,nrow=2,ncol=2,labels = letters[1:4]))
dev.off()




`9606.protein.links.v11.0` <- read.csv("./data/9606.protein.physical.links.full.v11.0.txt", sep="")
human_protein=as.data.frame(table(`9606.protein.links.v11.0`[,1]))
rownames(human_protein)<-substr(human_protein$Var1,6,20)

genes_infor$PPI_number_log=log10(human_protein[genes_infor$ensembl_peptide_id,2])
genes_infor$hg38_3UTR_length_log=log10(genes_infor$hg38_3UTR_length)

fimo <- read.delim("./data/meme_RNA_3/fimo.tsv", comment.char="#")
geneid=substr(fimo$sequence_name,17,31)
genes_infor$motif=rep("C-rich_motif-",length(genes_infor$ensembl_gene_id))
genes_infor$motif[genes_infor$ensembl_gene_id %in% geneid]="C-rich_motif+"
temp_x<-table(geneid)
genes_infor$motif_count=as.numeric(temp_x[genes_infor$ensembl_gene_id])
genes_infor$motif_count_density=log10(genes_infor$motif_count/genes_infor$hg38_3UTR_length)
genes_infor1=genes_infor[which(genes_infor$Gene.Type !="Oncogene"),]
GTEXmedian_tpm$motif_count_density<-genes_infor[GTEXmedian_tpm$Gene.ID,"motif_count_density"]
GTEXmedian_tpm$motif<-genes_infor[GTEXmedian_tpm$Gene.ID,"motif"]

p.ff_GTEx_length<-plotForest(GTEXmedian_tpm,"norm_exp","len_3UTR_log","Pearson's r(Log10 3'UTR Length VS\nGene Expression) in Normal samples(GTEx Project)",0.5,0.9,"1mm")

comp=list(c("C-rich_motif-","C-rich_motif+"))
ggplot(data = GTEXmedian_tpm,aes(y=norm_exp,x=motif,fill=motif))+geom_boxplot()+geom_signif(comparisons = comp,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1)+ylab("Log2(TPM+1) expression ")+facet_wrap(~cancer)
sig=c(rep("***",53))


comp=list(c("C-rich_motif-","C-rich_motif+"))
s1_motif=ggplot(data = GTEXmedian_tpm,aes(y=norm_exp,x=cancer,fill=motif))+geom_boxplot(aes(col=I("gray90")))+scale_fill_manual("Motif in 3'UTR", values=c("#76EE00","darkgreen"))+annotate("text",1:53,19,label=sig,col="red")+ylab("Gene Expression")+ggtitle("Normal samples")+ xlab("Dataset")+theme_hd()+theme(legend.position = "bottom",axis.text.x = element_text(angle = 90,hjust=0.95))



GTEXmedian_tpm$age_type2=genes_infor[GTEXmedian_tpm$Gene.ID,"age_type2"]
GTEX_cor_age_tsg<-get_cor_agetype(GTEXmedian_tpm,"norm_exp","len_3UTR_log","age_type2")

p_cor_age_gtex<-ggplot(data=GTEX_cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age class(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10 3'UTR Length VS\n Gene Expression)")+ggtitle("Normal samples\nGTEx Project RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))


p.ff_gtex_berbp<-plotForest(GTEXmedian_tpm,"norm_exp","rbp_site_binding_density2","Pearson's r(Log10(RBPs binding density \nin 3'UTR) VS Gene Expression) in Normal samples(GTEx Project)",0.5,0.9,"1mm")

GTEX_cor_age_tsg<-get_cor_agetype(GTEXmedian_tpm,"norm_exp","rbp_site_binding_density2","age_type2")

p_cor_age_gtex_rbp<-ggplot(data=GTEX_cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age class(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10(RBPs binding density \nin 3'UTR) VS Gene Expression)")+ggtitle("Normal samples\nGTEx Project RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))



GTEX_cor_age_tsg<-get_cor_agetype(GTEXmedian_tpm,"norm_exp","motif_count_density","age_type2")

p_cor_age_gtex_motif<-ggplot(data=GTEX_cor_age_tsg,aes(x=age,y=cor))+geom_boxplot(outlier.color = "white")+geom_jitter(width = 0.25,aes(col=Significance))+theme_hd()+xlab("Gene Age class(million years)")+theme(axis.text.x =element_text(  angle=30,hjust=1))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ylab("Pearson's r(Log10(RBPs binding density \nin 3'UTR) VS Gene Expression)")+ggtitle("Normal samples\nGTEx Project RNA-seq datasets")+scale_colour_manual(values = c("NS"="gray","*"="pink","**"="red","***"="darkred"))




p.ff_diffveri<-plotForest(all_data_other_cohort_filtered,"diff_exp","hg38_3UTR_length_log","Pearson's r(Log10 3'UTR Length VS abs(log2FC))\n Down-regulated TSG in Validation data",0.9,1.5,"1mm")

utr3m6a <- read.delim("./data/utr3m6a.txt", header=FALSE)
utr3m6a_MeTPeak=utr3m6a[utr3m6a[,8]=="MeTPeak",]
utr3m6a_MeTPeak$srx<-do.call(rbind,strsplit(utr3m6a_MeTPeak[,7],"_"))[,1]
repicm6a=unique(utr3m6a_MeTPeak[,c(6,7,11,12)])
repicm6a$gene=do.call(rbind,strsplit(repicm6a$V11,"\\."))[,1]
colnames(repicm6a)<-c("dataset","sample","genev","srx","gene")

all_srp <- read.delim("./data/all_srp.tsv")
all_srp_selcted=all_srp[all_srp$experiment_accession %in% utr3m6a_MeTPeak$srx,]

all_srp_selcted$cell=all_srp_selcted$cell.type
all_srp_selcted$cell[all_srp_selcted$cell==""]=all_srp_selcted$cell.line[all_srp_selcted$cell==""]
all_srp_selcted$cell[all_srp_selcted$cell==""]=all_srp_selcted$source_name[all_srp_selcted$cell==""]
all_srp_unique<-unique(all_srp_selcted[,c("experiment_accession","cell")])
rownames(all_srp_unique)<-all_srp_unique[,1]
repicm6a$cell=all_srp_unique[repicm6a$srx,"cell"]
repicm6a$cancer=paste(repicm6a$dataset,repicm6a$cell,sep="       ")
seleted_cell=c("embryonic kidney cells","embryonic stem (ES) cell","Endoderm","HEK293A-TOA"," Hek293T","HEK293T","HeLa","human embryonic stem cells","human neural progenitor cells","Jurkat","liver cell line","lymphoblastoid cell line (LCL)","OKMS inducible fibroblasts","organoid tissue","PCW11 human brain cortex","Primary CD4+ T","TREX")
repicm6a<-repicm6a[repicm6a$cell %in% seleted_cell,]
repicm6a_3utr<-c()
for(t1 in unique(repicm6a$cancer)){
  temp_genes_infor=genes_infor
  temp=repicm6a[repicm6a$cancer==t1,]
  x_count=table(temp$gene)
  select_genes<-names(x_count[x_count>1])
  if(length(select_genes)>1){
    temp_genes_infor$m6a=ifelse(temp_genes_infor$ensembl_gene_id %in% temp[,"gene"],"+","-")
    repicm6a_3utr<-rbind(repicm6a_3utr,cbind(cancer=t1,temp_genes_infor))
  }
  
}


tissue=c("Brain","Heart","Lung","Kidney","Liver","Muscle","Placenta","Stomach")
m6a_3utr<-c()
m6a_5utr<-c()
for(t1 in tissue){
  temp <- read.delim(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/result_UTR/m6a/GSE114150_",t1,"_3UTR.bed",sep=""), header=FALSE, stringsAsFactors=FALSE)
  temp_genes_infor=genes_infor
  temp_genes_infor$m6a=ifelse(temp_genes_infor$ensembl_gene_id %in% temp[,4],"+","-")
  m6a_3utr<-rbind(m6a_3utr,cbind(cancer=paste("GSE114150     ",t1),temp_genes_infor))
  
  temp2 <- read.delim(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/result_UTR/m6a/GSE114150_",t1,"_5UTR.bed",sep=""), header=FALSE, stringsAsFactors=FALSE)
  temp_genes_infor2=genes_infor
  temp_genes_infor2$m6a=ifelse(temp_genes_infor2$ensembl_gene_id %in% temp2[,4],"+","-")
  m6a_5utr<-rbind(m6a_5utr,cbind(cancer=paste("GSE114150       ",t1),temp_genes_infor2))
  
}

m6a_3utr2=m6a_3utr[m6a_3utr$Gene.Type!="Oncogene" ,]
m6a_5utr2=m6a_5utr[m6a_5utr$Gene.Type!="Oncogene" ,]

all_m6a<-rbind(m6a_3utr2,repicm6a_3utr[,colnames(m6a_3utr2)])


sif_repicm6a3_array=getsigm6a_diff(all_data_other_cohort_filtered,all_m6a,2,"diff_exp")
pff_array_repicm6a<-plotForesttwogroup(sif_repicm6a3_array,"Abs(Log2FC) of Down-regulated TSG in validation cohort","1mm","Potential m6A+","Potential m6A-",spacingx=1.4)






all_data_other_cohort_filtered$PPI_number_log<-genes_infor[all_data_other_cohort_filtered$ensembl_gene_id,"PPI_number_log"]



pff_ppi_veri<-plotForest(all_data_other_cohort_filtered,"diff_exp","PPI_number_log","Pearson's r(Log10 PPI number VS abs(log2FC))\n Down-regulated TSG in Validation data",0.9,1,"1mm")



x=table(all_m6a[all_m6a$m6a=="+",]$ensembl_gene_id)
all_data_other_cohort_filtered$m6a=0
all_data_other_cohort_filtered$m6a=x[all_data_other_cohort_filtered$ensembl_gene_id]
all_data_other_cohort_filtered$m6a[is.na(all_data_other_cohort_filtered$m6a)]=0
all_data_other_cohort_filtered$m6a_in<-ifelse(all_data_other_cohort_filtered$m6a>=2,"+","no")
all_data_other_cohort_filtered$m6a_in[all_data_other_cohort_filtered$m6a==0]<-"-"

p_m6adiff<-ggplot(data=all_data_other_cohort_filtered[which(all_data_other_cohort_filtered$cancer=="ccRCC_GSE53000_array" &all_data_other_cohort_filtered$m6a_in!="no"),],aes(x=m6a_in,y=diff_exp,fill=m6a_in))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("+","-")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold")+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+scale_x_discrete(labels= c("m6A-","m6A+"))+ggtitle("ccRCC_GSE53000_array")+ theme(plot.title = element_text(face = "bold"))+xlab("Potential m6A \n modification in 3'UTR")


pall<-ggarrange(p_cor_age_gtex,p_cor_age_gtex_rbp,p_cor_age_gtex_motif,p_m6adiff,labels = letters[c(3:5,8)], nrow=1,ncol=4,font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom",widths = c(1,1,1,0.8))


pall0<-ggarrange(p.ff_GTEx_length,p.ff_gtex_berbp,s1_rbp,labels = letters[1:3], nrow=1,ncol=3,font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1.3,1))

p_mga1<-ggarrange(p_cor_age_gtex,p.ff_diffveri,pff_array_repicm6a$p1,nrow=1,ncol=3,labels = letters[c(4:6)],font.label = list(size = 18, color = "black", face = "bold"),widths = c(1.2,1.2,2))
p_all1<-ggarrange( pall0,p_mga1,nrow=2,ncol=1,font.label = list(size = 18, color = "black", face = "bold"),heights  = c(2.2,1))


tiff("./result/figure_6_3UTR.tiff",width = 24.5,height = 17,res=300,units="in",compression = "lzw")
print(p_all1)
dev.off()


p_mga1<-ggarrange(pff_rbp_in,pff_ppi_veri,nrow=1,ncol=2,labels = letters[c(2:3)],font.label = list(size = 18, color = "black", face = "bold"),widths = c(1.37,1))
p_all1<-ggarrange( s1_motif,p_mga1,nrow=2,ncol=1,labels = c(letters[1],""),font.label = list(size = 18, color = "black", face = "bold"),heights =  c(2,1))


tiff("./result/supplymentary_figure13.tiff",width = 17.5,height = 13,res=300,units="in",compression = "lzw")
print(p_all1)
dev.off()



tiff("./result/supplymentary_figure14.tiff",width = 13,height = 13,res=300,units="in",compression = "lzw")
sp1=plot_scatter_2(all_data_other_cohort_filtered,"hg38_3UTR_length_log","diff_exp","Log10(Length of 3'UTR)",
                   "Abs(Log2FC) in Down-regulated TSGs",0.01,15,"black")+ggtitle("Validation data")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4,scales = "free_y")
print(sp1)

dev.off()


tiff("./result/supplymentary_figure15.tiff",width = 13,height = 13,res=300,units="in",compression = "lzw")
sp1=plot_scatter_2(all_data_other_cohort_filtered,"PPI_number_log","diff_exp","Log10(PPI number)",
                   "Abs(Log2FC) in Down-regulated TSGs",0.01,15,"black")+ggtitle("Validation data")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4,scales = "free_y")
print(sp1)

dev.off()


