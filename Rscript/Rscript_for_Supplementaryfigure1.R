

source('functions/functions_evolution.r')
load("./data/genes_infor_downloaded_from_ensembl97.rdata")
load("./data/other_species_utr_more.rdata")
gene_age_infor <- read.delim("./data/gene_age_infor.txt")
rownames(gene_age_infor)<-gene_age_infor$Age.class

trans_infor<-unique(all_data_genes$get_gene_pos[,c("ensembl_transcript_id","transcript_appris")])
rownames(trans_infor)<-trans_infor[,1]



get_diff_sig<-function(all_data,col.x){
  sp=unique(all_data$species)
  p_value=c()
  for(i in sp){
    temp.data<-all_data[all_data$species==i,]
    x <- try(wilcox.test(temp.data[,col.x]~temp.data$Gene.Type),TRUE)
    if (class(x) == "try-error") {
      p_value=c(p_value,NA)
    }else{
      p_value=c(p_value,x$p.value)
    }
  }
  names(p_value)<-sp
  return(p_value)
}

plotForSpecies<-function(all_data,col.x,lab.x){
  wilcox_p<-get_diff_sig(all_data,col.x)
  wilcox_Sig <- ifelse(wilcox_p>0.05,"NS",ifelse(wilcox_p >0.01, "*" , 
                                                 ifelse(wilcox_p >0.001, "**" , "***")))
  max.x<-quantile(all_data[,col.x],1,na.rm = T)
  
  p=ggplot(data=all_data[!is.na(all_data$Gene.Type),],aes(x=species,y=get(col.x),fill=Gene.Type))+geom_boxplot(width=0.4,aes(col=I("gray80")))+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+xlab("  ")+annotate("text",x=1:34,y=rep(max.x,34),label=wilcox_Sig[levels(all_data$species)],col="red",fontface =2)+theme(axis.text.x = element_text(angle = 90,hjust = 0.95))
  return(p+ylab(lab.x))
  #+stat_summary(fun.data =  countFunction,aes(group=Gene.Type), geom="text", size = 3,col="black", hjust = 0.5, position = position_dodge(0.6),fontface =2)
}


path <- "./data/Species_we_selected_V2.csv"
Data<-read.csv(path)
Species.select<- Data$Species
gene_pos_data_human=cbind(genes_infor[which(genes_infor$Gene.Type!="Oncogene"),c("ensembl_gene_id","genes_trans_length","hg38_5UTR_length","hg38_3UTR_length","Gene.Type")],species="Human")
colnames(gene_pos_data_human)<-colnames(result.utr.others)
all_data=rbind(gene_pos_data_human,result.utr.others)
all_data$len_3UTR[all_data$len_3UTR==0]=NA
all_data$len_5UTR[all_data$len_5UTR==0]=NA
level_x=Species.select[34:1]
#all_data$species<-factor(all_data$species,levels = level_x[7:1])

all_data$hg38_3UTR_length_log<-log10(all_data$len_3UTR)
all_data$hg38_5UTR_length_log<-log10(all_data$len_5UTR)
all_data$cds_length_log<-log10(all_data$trans_length-all_data$len_5UTR-all_data$len_3UTR)
all_data$hg38_3utr_ratio_all<-all_data$len_3UTR/all_data$trans_length
all_data$hg38_3utr_ratio_non5UTR<-all_data$len_3UTR/(all_data$trans_length-all_data$len_5UTR)
all_data$hg38_3utr_to_cds<-all_data$hg38_3UTR_length_log-all_data$cds_length_log
all_data$hg38_5utr_to_cds<-all_data$hg38_5UTR_length_log-all_data$cds_length_log
all_data$log_ratio_3UTR_5UTR<-log10(all_data$len_3UTR/all_data$len_5UTR)
all_data$hg38_5utr_ratio_all<-all_data$len_5UTR/all_data$trans_length
all_data$genes_trans_length_log<-log10(all_data$trans_length)

all_data$species<-factor(all_data$species,levels = level_x)
s0<-plotForSpecies(all_data,"hg38_3UTR_length_log","Log10(Length of 3'UTR)")
s1<-plotForSpecies(all_data,"hg38_3utr_ratio_all","Normalized length of 3'UTR")
s0= s0+theme(axis.text.x = element_blank())+ scale_y_continuous(breaks=c(1.0,2.0, 3.0, 4.0),labels = c("1.00","2.00", "3.00", "4.00"))
###################################################################33

genes_infor$transcript_appris=trans_infor[genes_infor$genes_main_transcript,2]
genes_infor_2=genes_infor[genes_infor$Gene.Type!="Oncogene",]
gene_paralog=all_data_genes$gene_paralog
gene_paralog=gene_paralog[gene_paralog$hsapiens_paralog_ensembl_gene!="",]
TSG_genes=genes_infor_2[genes_infor_2$Gene.Type=="TSG",]$ensembl_gene_id
TSG_paralog_genes=gene_paralog[gene_paralog$ensembl_gene_id %in% TSG_genes & 
                                 gene_paralog$hsapiens_paralog_ensembl_gene %in% genes_infor_2$ensembl_gene_id[genes_infor_2$Gene.Type=="Non-Cancer"],]
TSG_paralog_non_cancer=genes_infor_2[unique(as.character(c(TSG_genes,TSG_paralog_genes$hsapiens_paralog_ensembl_gene))),]  

TSG_paralog_non_cancer$Gene.Type[TSG_paralog_non_cancer$Gene.Type=="Non-Cancer"]="Non_Cancer_paralogs"

p1_3utr_prin<-plot_non_paired_paralog(TSG_paralog_non_cancer[ TSG_paralog_non_cancer$transcript_appris=="principal1",],"hg38_3UTR_length_log","Log10(Length of 3'UTR)")+ggtitle("APPRIS-principal1\n paralogs")


p1_3utr_prin0=plot_TSG_non_cancer(genes_infor[ genes_infor$transcript_appris=="principal1" &genes_infor$Gene.Type!="Oncogene",],"hg38_3UTR_length_log","Log10(Length of 3'UTR)")+ggtitle("APPRIS-principal1")


p1_3utr_nonprin<-plot_non_paired_paralog(TSG_paralog_non_cancer[ TSG_paralog_non_cancer$transcript_appris!="principal1",],"hg38_3UTR_length_log","Log10(Length of 3'UTR)")+ggtitle("Non-APPRIS-principal1\n paralogs")

p1_3utr_nonprin0=plot_TSG_non_cancer(genes_infor[ genes_infor$transcript_appris!="principal1" &genes_infor$Gene.Type!="Oncogene",],"hg38_3UTR_length_log","Log10(Length of 3'UTR)")+ggtitle("Non-APPRIS-principal1")

cut_x=c(0,7,14,18,23,27)
cut_lable=c()
for(i in 2:length(cut_x)){
  x=cut_x[i-1]
  if(x==0){x=1}
  cut_lable<-c(cut_lable,paste("(",gene_age_infor[x,4]," - ",gene_age_infor[cut_x[i],4],"]",sep=""))
}
cut_lable[5]=">1119.25 "
p_class<-plot_classified_hist(genes_infor$Age_class[!is.na(genes_infor$Age_class)],cut_x+0.5,cut_lable,"Age class")+theme_hd()
genes_infor<-genes_infor %>% mutate(age_type2=cut(genes_infor$Age_class, breaks=cut_x,labels = cut_lable))
rownames(genes_infor)<-genes_infor$ensembl_gene_id

cols_x=brewer.pal(n = 8, name = "Set2")
p_class_num<-ggplot(data=genes_infor[!is.na(genes_infor$age_type2),],aes(x=age_type2,fill=age_type2))+geom_bar()+theme_hd()+xlab("Gene age group")+theme(axis.text.x =element_text(  angle=90,hjust=1))+ylab("Number of genes")+scale_fill_manual(name="Class",values =   cols_x[1:5])+guides(fill=guide_legend(nrow=2,byrow=TRUE))

p_class_all<-ggarrange(p_class,p_class_num,nrow=1,ncol=2, font.label = list(size = 18, color = "black", face = "bold"),labels = letters[5:6],common.legend = T,legend = "bottom")

p_princ<-ggarrange(p1_3utr_prin0+ylim(0,5),p1_3utr_nonprin0+ylim(0,5),p1_3utr_prin+ylim(0,5),p1_3utr_nonprin+ylim(0,5),nrow=2,ncol=2, font.label = list(size = 18, color = "black", face = "bold"),labels = letters[1:4])



p_all0<-ggarrange(p_princ,p_class_all,nrow=1,ncol=2, font.label = list(size = 18, color = "black", face = "bold"))


palls1=ggarrange(p_all0,s0,s1,nrow=3,ncol=1,labels = c("",letters[7:8]), font.label = list(size = 18, color = "black", face = "bold"),heights  = c(1.4,1,1.3),common.legend = T,legend = "bottom")

tiff("./result/supplymentary_figure1.tiff",width = 13,height = 14,res=300,units="in",compression = "lzw")
print(palls1)
dev.off()




