setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/final_program_UTR9_21")
load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')



`9606.protein.links.v11.0` <- read.csv("./data/9606.protein.physical.links.full.v11.0.txt", sep="")
human_protein=as.data.frame(table(`9606.protein.links.v11.0`[,1]))
rownames(human_protein)<-substr(human_protein$Var1,6,20)

genes_infor$PPI_number_log=log10(human_protein[genes_infor$ensembl_peptide_id,2])
genes_infor$hg38_3UTR_length_log=log10(genes_infor$hg38_3UTR_length)
genes_infor$rbp<-ifelse(is.na(genes_infor$rbp_site_binding_density2),"RBP-","RBP+")
all_data_tsg_filtered$rbp=genes_infor[all_data_tsg_filtered$ensembl_gene_id,"rbp"]
all_data_tsg_filtered$motif=genes_infor[all_data_tsg_filtered$ensembl_gene_id,"motif"]

p.ff_rbp<-plotForest(all_data_tsg_filtered,"diff_exp","rbp_site_binding_density2","Pearson's r(Log10 RBP binding density VS abs(log2FC))\n Down-regulated TSG in TCGA",0.9,0.9,"1cm")


p.ff_diff<-plotForest(all_data_tsg_filtered,"diff_exp","hg38_3UTR_length_log","Pearson's r(Log10_3UTR_Length VS abs(log2FC))\n Down-regulated TSG in TCGA",0.9,1,"5cm")

cor_plot_len_diff<-get_cor_for_each_cancer_all(all_data_tsg_filtered,"diff_exp","hg38_3UTR_length_log","pearson")

cancer="COAD"
cancer_temp=all_data_tsg_filtered[all_data_tsg_filtered$cancer=="COAD",]
cancer_temp$diff_exp=abs(cancer_temp$diff_exp)
sd1=plot_scatter_2(cancer_temp,"hg38_3UTR_length_log","diff_exp","Log10(Length of 3'UTR)",
                   "Abs(Log2FC) of Down-regulated TSG",0,4.5,"black")+ggtitle(paste("TCGA-COAD"))+theme(plot.title = element_text(hjust = 0.5)) +theme_hd()
print(sd1)



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
seleted_cell=c("embryonic kidney cells","embryonic stem (ES) cell","Endoderm","HEK293A-TOA"," Hek293T","HEK293T","human embryonic stem cells","human neural progenitor cells","Jurkat","liver cell line","lymphoblastoid cell line (LCL)","OKMS inducible fibroblasts","organoid tissue","PCW11 human brain cortex","Primary CD4+ T","TREX")
repicm6a<-repicm6a[repicm6a$cell %in% seleted_cell,]
repicm6a_3utr<-c()

for(t1 in unique(repicm6a$cancer)){
  temp_genes_infor=genes_infor
  temp=repicm6a[repicm6a$cancer==t1,]
  x2=unique(temp[,1:2])
  x_count=table(temp$gene)
  select_genes<-names(x_count[x_count>1])
  if(length(select_genes)>1){
    temp_genes_infor$m6a=ifelse(temp_genes_infor$ensembl_gene_id %in% temp[,"gene"],"+","-")
    repicm6a_3utr<-rbind(repicm6a_3utr,cbind(cancer=t1,temp_genes_infor,count=dim(x2)[1]))
  }
  
}

sig_rbp<-get_diff(all_data_tsg_filtered,"rbp","diff_exp","RBP+")
pff_rbp_all<-plotForesttwogroup(sig_rbp,"Abs(Log2FC) of Down-regulated TSG in TCGA","4.2mm","RBP+","RBP-")
pff_rbp_in<-pff_rbp_all$p1

p_rbp_box<-ggplot(data=all_data_tsg_filtered[all_data_tsg_filtered$cancer=="COAD",],aes(x=rbp,y=diff_exp,fill=rbp))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("RBP+","RBP-")), test='wilcox.test', map_signif_level=T,col="black",textsize = 4.8,fontface="bold",y_position = 4.3)+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c("#76EE00","darkgreen"))+theme_hd()+ theme(plot.title = element_text(face = "bold"))+xlab("RBP binding site in 3'UTR")+ggtitle("TCGA-COAD")



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

sif_repicm6a3=getsigm6a_diff(all_data_tsg_filtered,all_m6a,2,"diff_exp")
pff_tcga_repicm6a_all<-plotForesttwogroup(sif_repicm6a3,"Abs(Log2FC) of Down-regulated TSG in TCGA","4.2mm","Potential m6A+","Potential m6A-")
pff_tcga_repicm6a<-pff_tcga_repicm6a_all$p1
ggarrange(pff_tcga_repicm6a)

sif_repicm6a3_len=getsigm6a_difflen(all_m6a[all_m6a$Gene.Type=="TSG",],4,"hg38_3UTR_length_log","m6a")
pff_repicm6a_len_all<-plotForesttwogroup(sif_repicm6a3_len,"Log10(3UTR Length) of TSGs","1mm","3UTR m6A+","3UTR m6A-")
pff_repicm6a_len<-pff_repicm6a_len_all$p1
temp_Example<-repicm6a_3utr[repicm6a_3utr$Gene.Type=="TSG" & repicm6a_3utr$cancer=="SRP007335       embryonic kidney cells",]
rownames(temp_Example)<-temp_Example$ensembl_gene_id
p_m6alen<-ggplot(data=temp_Example,aes(x=m6a,y=hg38_3UTR_length_log,fill=m6a))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("+","-")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold",y_position = 4.4)+ylab("Log10(3UTR Length) of TSGs") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+scale_x_discrete(labels= c("m6A-","m6A+"))+ggtitle("SRP007335 \nembryonic kidney cells")+ theme(plot.title = element_text(face = "bold"))+xlab("m6A modification in 3'UTR")


x=table(all_m6a[all_m6a$m6a=="+",]$ensembl_gene_id)
all_data_tsg_filtered$m6a=0
all_data_tsg_filtered$m6a=x[all_data_tsg_filtered$ensembl_gene_id]
all_data_tsg_filtered$m6a[is.na(all_data_tsg_filtered$m6a)]=0
all_data_tsg_filtered$m6a_in<-ifelse(all_data_tsg_filtered$m6a>=2,"+","no")
all_data_tsg_filtered$m6a_in[all_data_tsg_filtered$m6a==0]<-"-"

p_m6adiff<-ggplot(data=all_data_tsg_filtered[which(all_data_tsg_filtered$cancer=="COAD" &all_data_tsg_filtered$m6a_in!="no"),],aes(x=m6a_in,y=diff_exp,fill=m6a_in))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("+","-")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold",y_position = 4.2)+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+scale_x_discrete(labels= c("m6A-","m6A+"))+ggtitle("TCGA-COAD")+ theme(plot.title = element_text(face = "bold"))+xlab("Potential m6A modification in 3'UTR")

all_data_tsg_filtered$rbp_site_binding2=log10(genes_infor[all_data_tsg_filtered$ensembl_gene_id,"rbp_site_binding2"])
all_data_tsg_filtered$rbp_type_count=log10(genes_infor[all_data_tsg_filtered$ensembl_gene_id,"rbp_type_count"])

pcor1=plotForest(all_data_tsg_filtered,"diff_exp","rbp_site_binding2","Pearson's r(Log10 No.of RBP binding site VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")
pcor2=plotForest(all_data_tsg_filtered,"diff_exp","rbp_type_count","Pearson's r(Log10 No.of RBPs VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")
pcor3=plotForest(all_data_tsg_filtered,"diff_exp","rbp_site_binding_density2","Pearson's r(Log10 rbp_site_binding_density VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1.1,"1mm")



p1=ggplot(data=all_data_tsg_filtered,aes(x=rbp,y=diff_exp,fill=rbp))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("RBP+","RBP-")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold")+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+ theme(plot.title = element_text(face = "bold"),axis.text.x = element_text(angle = 90))+xlab("RBP binding site in 3'UTR")+facet_wrap(~cancer,nrow=1)


p2=ggplot(data=all_data_tsg_filtered,aes(x=motif,y=diff_exp,fill=motif))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+geom_signif(comparisons = list(c("With","Without")), test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold")+ylab("Abs(Log2FC) of Down-regulated TSG") +scale_fill_manual(values=c( "#F0FFFF","#1874CD"))+theme_hd()+ theme(plot.title = element_text(face = "bold"),axis.text.x = element_text(angle = 90))+xlab("C-rich motif in 3'UTR")+facet_wrap(~cancer,nrow=1)



tiff("./result/test_figure_4_3UTR.tiff",width = 20,height = 12,res=300,units="in",compression = "lzw")
print(ggarrange(pcor1,pcor2,pcor3,p1,nrow=2,ncol=2,labels = letters[1:4]))
dev.off()



cut_tcga<-getSMDforestcutoff(all_data_tsg_filtered,all_m6a,"diff_exp","Potential  m6A+","Potential  m6A-")

ggplot(data = cut_tcga,aes(x=cutoff,y=smd,col=cancer))+geom_boxplot()


all_data_tsg_filtered$PPI_number_log<-genes_infor[all_data_tsg_filtered$ensembl_gene_id,"PPI_number_log"]

ppi_TSG<-plot_scatter_2(genes_infor[which(genes_infor$Gene.Type=="TSG"),],"hg38_3UTR_length_log","PPI_number_log","Log10(Length of 3'UTR)","Log10 PPI number",1,20,"black")+ggtitle("TSG")+ylim(1.3,3.8)+theme_hd()


pff_PPI<- plotForest(all_data_tsg_filtered,"diff_exp","PPI_number_log","Pearson's r(Log10 PPI number VS abs(log2FC))\n Down-regulated TSG in TCGA datasets",0.9,1,"1mm")

cor_plot_ppi_diff<-get_cor_for_each_cancer_all(all_data_tsg_filtered,"diff_exp","PPI_number_log","pearson")
cor_plot_rbp_diff<-get_cor_for_each_cancer_all(all_data_tsg_filtered,"diff_exp","rbp_site_binding_density2","pearson")
cor_plot_rbp_count_diff<-get_cor_for_each_cancer_all(all_data_tsg_filtered,"diff_exp","rbp_site_binding2","pearson")

pall0<-ggarrange(p.ff_diff,sd1,labels = letters[1:2], nrow=1,ncol=2,font.label = list(size = 18, color = "black", face = "bold"),widths = c(3,2))
pall1<-ggarrange(pff_rbp_in,p_rbp_box,p_m6alen,labels = letters[3:5], nrow=1,ncol=3,font.label = list(size = 18, color = "black", face = "bold"),widths = c(3,1,1))

p_mga0<-ggarrange(pff_repicm6a_len,p_m6adiff,ppi_TSG,nrow=1,ncol=3,labels = c(letters[c(6:7,9)]),font.label = list(size = 18, color = "black", face = "bold"),widths = c(3,1,1))
p_mga1<-ggarrange(pff_tcga_repicm6a,pff_PPI,nrow=1,ncol=2,labels = letters[c(8,10)],font.label = list(size = 18, color = "black", face = "bold"),widths = c(3,2))

p_mga2<-ggarrange(pall0,pall1,p_mga0,p_mga1, nrow=4,ncol=1, font.label = list(size = 18, color = "black", face = "bold"))

tiff("./result/figure_4_3UTR.tiff",width = 19,height = 19,res=300,units="in",compression = "lzw")
print(p_mga2)
dev.off()



library(readxl )
percentage_cancer <-as.data.frame( read_excel("./data/percentage_Comprehensive Characterization of Cancer Driver Genes and Mutations.xlsx"))
percentage_cancer$cancer[percentage_cancer$cancer=="COADREAD"]="COAD"
percentage_cancer=percentage_cancer[-which(percentage_cancer$cancer=="STAD"),]
rownames(percentage_cancer)<-percentage_cancer$cancer
percentage_cancer$type="Unknown"
percentage_cancer$type[as.numeric(percentage_cancer$TSG)>=0.60]="TSG driver cancers"
percentage_cancer$type[as.numeric(percentage_cancer$Oncogene)>=0.60]="Oncogene driver cancers"
percentage_cancer$cancer=factor(percentage_cancer$cancer,levels = percentage_cancer$cancer[order(percentage_cancer$TSG)])
pd0<-ggplot(data=percentage_cancer,aes(x=cancer,y=TSG,fill=type))+geom_bar(stat="identity",position = "dodge")+coord_flip()
all_data_tsg_filtered$FC.larger.than2<- ifelse(all_data_tsg_filtered$diff_exp>1,"FC > 1 ", "FC < 1")
x<-table(all_data_tsg_filtered$FC.larger.than2,all_data_tsg_filtered$cancer)
ratio.x<-x[2,]/apply(x,2,sum)
ratio.x<-sort(ratio.x)
data.x<-data.frame(cancer=names(ratio.x),ratio=ratio.x,n=x[2,names(ratio.x)])
data.x$cancer<-factor(data.x$cancer,levels=names(ratio.x))
data.x$percent_cancer_driver<-percentage_cancer[as.character(data.x$cancer),"TSG"]
rownames(data.x)<-as.character(data.x$cancer)



all_cor=rbind(cbind(type="3utr_diff",cor_plot_len_diff),cbind(type="ppi_diff",cor_plot_ppi_diff),cbind(type="rbp_diff",cor_plot_rbp_diff),cbind(type="rbp_count_diff",cor_plot_rbp_count_diff))

all_cor$percent_cancer_driver<-percentage_cancer[as.character( all_cor$Cancer),2]
all_cor$ratio_higher_TSG<-data.x[as.character( all_cor$Cancer),"ratio"]
all_cor$number_higher_TSG<-data.x[as.character( all_cor$Cancer),"n"]
all_cor$cancer_type<-percentage_cancer[as.character( all_cor$Cancer),"type"]
all_cor$r2<-all_cor$cor^2

pr2=ggscatter(all_cor[all_cor$type=="ppi_diff",] ,  y = "r2", x = "ratio_higher_TSG",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's R2(Log10 PPI number\nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")


pd1 <- ggscatter(all_cor[all_cor$type=="3utr_diff",] ,  y = "cor", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's r(Log10 3'UTR Length\nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")
#+ylim(0,0.1)+xlim(0.15,0.9)

pd2_count <- ggscatter(all_cor[all_cor$type=="rbp_count_diff",] ,  y = "r2", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's R2(Log10 No.of RBP binding site\nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")+ylim(-0.02,0.11)+xlim(0.15,0.9)

pd2_denity <- ggscatter(all_cor[all_cor$type=="rbp_diff",] ,  y = "r2", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's R2(Log10 RBP binding density \nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")+ylim(-0.02,0.11)+xlim(0.15,0.9)

pd2 <- ggscatter(all_cor[all_cor$type=="ppi_diff",] ,  y = "cor", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "Cancer")+xlab("Percentage of TSG driver genes")+ylab("Pearson's r(Log10 PPI number\nVS abs(log2FC))")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",label.y = 0.06,aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")



cut_off_tcga<-as.data.frame(cut_tcga[cut_tcga$cutoff==">= 2" ,])
cut_off_tcga$smd<-abs(as.numeric(cut_off_tcga$smd))
cut_off_tcga$percent_cancer_driver<-percentage_cancer[as.character(cut_off_tcga$cancer),2]
cut_off_tcga$ratio_higher_TSG<-data.x[as.character( cut_off_tcga$cancer),"ratio"]
cut_off_tcga$number_higher_TSG<-data.x[as.character( cut_off_tcga$cancer),"n"]
cut_off_tcga$cancer_type<-percentage_cancer[as.character( cut_off_tcga$cancer),"type"]

smd_cor<-data.frame(data.x,m6a=pff_tcga_repicm6a_all$SMD[as.character(data.x$cancer)],rbp=pff_rbp_all$SMD[as.character(data.x$cancer)])

pd_m6a<-ggscatter(smd_cor ,  y = "m6a", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "cancer")+xlab("Percentage of TSG driver genes")+ylab("SMD of Abs(log2FC) \nbewtween m6A+ vs m6A-")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")

pd_rbp<-ggscatter(smd_cor ,  y = "rbp", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "cancer")+xlab("Percentage of TSG driver genes")+ylab("SMD of Abs(log2FC) \nbewtween RBP+ vs RBP-")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")


pd3 <- ggscatter(cut_off_tcga ,  y = "smd", x = "percent_cancer_driver",alpha=I(0.8),add = "reg.line",  label = "cancer")+xlab("Percentage of TSG driver genes")+ylab("Abs(SMD of Abs(log2FC) \nbewtween m6A+ vs m6A-)")+ theme_hd()+stat_cor(method = "spearman",fontface="bold",output.type="text",aes(label=paste(sub("R = ","Rho=",..r.label..),", ",..p.label..,sep="")))+ggtitle("Down-regulated TSG")+ylim(0.2,1.3)+xlim(0.15,0.9)



pd<-ggarrange(pd1,pd2,pd_rbp,pd_m6a,nrow=2,ncol=2,labels = letters[1:4],font.label = list(size = 18, color = "black", face = "bold"))


tiff("./result/figure_5_3UTR.tiff",width = 16,height = 16,res=300,units="in",compression = "lzw")
print(pd)
dev.off()


sample_num_count=as.data.frame(sample_num_all)
colnames(sample_num_count)<-c("cancer","type","sample_id")
exp_sample_sum=as.data.frame(sort(table(sample_num_all[sample_num_all[,2]=="exp",1])))

p1<-ggplot(data=exp_sample_sum, aes(x=Var1, y=Freq,fill=I("#00BFFF"))) +
  geom_bar(stat="identity")+ geom_text(aes(label = Freq, y = Freq))+xlab("TCGA Datasets")+ylab("Number of paired tumor and\n tumor adjacent normal samples")+theme_hd()+coord_flip()
temp_count<-unique(repicm6a_3utr[,c("cancer","count")])
p2<-ggplot(data=temp_count[order(temp_count$count),], aes(x=cancer, y=count,fill=I("#00BFFF"))) +geom_bar(stat="identity")+ geom_text(aes(label = count, y = count))+xlab("REPIC m6A Datasets")+ylab("Number of samples")+theme_hd()+coord_flip()+theme(axis.text.y = element_text(hjust=0))+ylim(0,13)
tiff("./result/supplymentary_figure2.tiff",width = 10,height = 10,res=300,units="in",compression = "lzw")

p=ggarrange(p1,p2,nrow=2,ncol=1,labels = letters[1:2],font.label = list(size = 18, color = "black", face = "bold"))
print(p)
dev.off()


#supplymentary_


p_cor_rbp_tsg<-plotForest(all_data_tsg_filtered,"diff_exp","rbp_site_binding_density2","Pearson's r(Log10(RBPs binding density in 3'UTR) VS abs(log2FC))\n Down-regulated TSG in Validation data",0.9,1.1,"1cm")


tiff("./result/supplymentary_figure11.tiff",width = 13,height = 13,res=300,units="in",compression = "lzw")
sp1=plot_scatter_2(all_data_tsg_filtered[ all_data_tsg_filtered$cancer !="COAD",],"hg38_3UTR_length_log","diff_exp","Log10(Length of 3'UTR)",
                   "Abs(Log2FC) of Down-regulated TSGs",0.01,15,"black")+ggtitle("TCGA RNA-seq data")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4)
print(sp1)

dev.off()




tiff("./result/supplymentary_figure12.tiff",width = 13,height = 13,res=300,units="in",compression = "lzw")
sp1=plot_scatter_2(all_data_tsg_filtered[ all_data_tsg_filtered$cancer !="COAD",],"PPI_number_log","diff_exp","Log10(PPI number)",
                   "Abs(Log2FC) in Down-regulated TSGs",0.01,15,"black")+ggtitle("TCGA RNA-seq data")+theme(strip.background = element_blank())+ guides(colour = "none",fill = "none")+scale_color_manual(values=c("black"))+facet_wrap(~cancer,nrow=4)
print(sp1)

dev.off()






