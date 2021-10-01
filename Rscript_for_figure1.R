load("./data/combined_gene_infor_expression_data.RData")
source('functions/functions_evolution.r')
trans_infor<-unique(all_data_genes$get_gene_pos[,c("ensembl_transcript_id","transcript_appris")])
rownames(trans_infor)<-trans_infor[,1]

genes_infor_out_put<-genes_infor[,c("ensembl_gene_id","hgnc_symbol","genes_main_transcript","Gene.Type","hg38_3UTR_length","Age_class","Origin.time..million.years.ago..","age_type2")]

colnames(genes_infor_out_put)<-c("Ensembl_gebe_ID","HGNC_Symbol","Candidate_transcript","Gene_Type","3UTR length","Age Class","Original_time(Million years age","Age Group")

write.csv(genes_infor_out_put,file="./result/genes_infor.csv",row.names = F)





genes_infor$transcript_appris=trans_infor[genes_infor$genes_main_transcript,2]
genes_infor_2=genes_infor[genes_infor!="Oncogene",]
gene_paralog=all_data_genes$gene_paralog
gene_paralog=gene_paralog[gene_paralog$hsapiens_paralog_ensembl_gene!="",]
TSG_genes=genes_infor_2[genes_infor_2$Gene.Type=="TSG",]$ensembl_gene_id
TSG_paralog_genes=gene_paralog[gene_paralog$ensembl_gene_id %in% TSG_genes & 
                                 gene_paralog$hsapiens_paralog_ensembl_gene %in% genes_infor_2$ensembl_gene_id[genes_infor_2$Gene.Type=="Non-Cancer"],]
TSG_paralog_non_cancer=genes_infor_2[unique(as.character(c(TSG_genes,TSG_paralog_genes$hsapiens_paralog_ensembl_gene))),]  

TSG_paralog_non_cancer$Gene.Type[TSG_paralog_non_cancer$Gene.Type=="Non-Cancer"]="Non_Cancer_paralogs"



temp_gene_pro<-unique(all_data_genes$get_gene_pos[!is.na(all_data_genes$get_gene_pos$`3_utr_start`) & all_data_genes$get_gene_pos$transcript_biotype=="protein_coding",c("ensembl_gene_id","ensembl_transcript_id","3_utr_start","3_utr_end","len_3UTR")])

temp_gene_pro2<-unique(all_data_genes$get_gene_pos[!is.na(all_data_genes$get_gene_pos$`3_utr_start`) & all_data_genes$get_gene_pos$transcript_biotype=="protein_coding",c("ensembl_gene_id","ensembl_transcript_id","len_3UTR")])

temp_gene_pro3=unique(temp_gene_pro2[,c("ensembl_gene_id","len_3UTR")])
utr3_count<-table(temp_gene_pro3$ensembl_gene_id)
genes_infor$count_3UTR<-as.numeric(utr3_count[genes_infor$ensembl_gene_id])


temp_aver_3UTR<-aggregate( temp_gene_pro2$len_3UTR,list( temp_gene_pro2$ensembl_gene_id),mean)
rownames(temp_aver_3UTR)<-temp_aver_3UTR[,1]
genes_infor$aver_3UTR<-log10(temp_aver_3UTR[genes_infor$ensembl_gene_id,2])

temp_max_3UTR<-aggregate( temp_gene_pro2$len_3UTR,list( temp_gene_pro2$ensembl_gene_id),max)
rownames(temp_max_3UTR)<-temp_max_3UTR[,1]
genes_infor$max_3UTR<-log10(temp_max_3UTR[genes_infor$ensembl_gene_id,2])

temp_min_3UTR<-aggregate( temp_gene_pro2$len_3UTR,list( temp_gene_pro2$ensembl_gene_id),min)
rownames(temp_min_3UTR)<-temp_min_3UTR[,1]
genes_infor$min_3UTR<-log10(temp_min_3UTR[genes_infor$ensembl_gene_id,2])

gene_trans<-unique(all_data_genes$get_gene_pos[,1:2])
trans_count<-table(gene_trans[,1])
genes_infor$all_trans_count<-as.numeric(trans_count[genes_infor$ensembl_gene_id])

gene_trans_pro<-unique(all_data_genes$get_gene_pos[all_data_genes$get_gene_pos$transcript_biotype=="protein_coding",1:2])





trans_pro_count<-table(gene_trans_pro[,1])
genes_infor$all_trans_pro_count<-as.numeric(trans_pro_count[genes_infor$ensembl_gene_id])



p_trans_count=plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],"all_trans_count","No. of transcripts per gene")+scale_y_log10()
p_trans_pro_count=plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],"all_trans_pro_count","No. of transcripts per gene")+scale_y_log10()
p_utr3_count=plot_TSG_non_cancer(genes_infor[genes_infor$Gene.Type!="Oncogene",],"count_3UTR","No. of alternative 3'UTR per gene")+scale_y_log10()

p1_3utr<-plot_non_paired_paralog(TSG_paralog_non_cancer,"hg38_3UTR_length_log","Log10(Length of 3'UTR)")
p1_trans<-plot_non_paired_paralog(TSG_paralog_non_cancer,"genes_trans_length_log","Log10(Length of transcript)")
p1_3utr_all<-plot_non_paired_paralog(TSG_paralog_non_cancer,"hg38_3utr_ratio_all","3UTR Ratio")
genes_infor=genes_infor[genes_infor$Gene.Type!="Oncogene",]
p0_3utr<-plot_TSG_non_cancer(genes_infor,"hg38_3UTR_length_log","Log10(Length of 3'UTR)")
p0_trans<-plot_TSG_non_cancer(genes_infor,"genes_trans_length_log","Log10(Length of transcript)")
p0_3utr_all<-plot_TSG_non_cancer(genes_infor,"hg38_3utr_ratio_all","3'UTR Ratio")


p_aver_3UTR<-plot_TSG_non_cancer(genes_infor,"aver_3UTR","Log10(average 3'UTR length)")
p_max_3UTR<-plot_TSG_non_cancer(genes_infor,"max_3UTR","Log10(max 3'UTR length)")
p_min_3UTR<-plot_TSG_non_cancer(genes_infor,"min_3UTR","Log10(min 3'UTR length)")


genes_infor$age_type3<-as.character(genes_infor$age_type2)
genes_infor$age_type3[genes_infor$age_type3%in% c("(67.6 - 355.7]","(0 - 67.6]")]="(0 - 355.7]"
p1_3UTR_age<-plotAgelength(genes_infor,"hg38_3UTR_length_log","Log10(Length of 3'UTR)",0,"age_type3")

countFunction <- function(x){
  return(data.frame(y=-0.5,label=round(length(x),2)))}
p_age_3utr1=ggplot(genes_infor[!is.na(genes_infor$age_type3),],aes(x=age_type3,y=hg38_3UTR_length_log,fill=age_type3))+geom_boxplot(show.legend = FALSE)+geom_signif(comparisons=getComp(unique(genes_infor$age_type3[!is.na(genes_infor$age_type3)])), test='wilcox.test', map_signif_level=TRUE,col="black",step_increase=0.1) +theme_hd()+scale_fill_manual(values=brewer.pal(n = 5, name = 'Blues')[2:5])+theme(axis.text.x = element_text(angle = 90,hjust = 1))+ylab("Log10(3'UTR Length)")+xlab("Gene Age class\n(million years)")#+stat_summary(fun.data =  countFunction, geom="text", size = 4,col="black")

p_age_3utr0=plotSubSigOneSide(genes_infor[!is.na(genes_infor$age_type3),],"age_type3", "hg38_3UTR_length_log","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"Gene Age class\n(million years)","Log10(Length of 3'UTR)",side="less") + labs(fill = "Gene.Type")+theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")




#####################################
cancermine_collated <- read.delim("./data/cancermine_collated.tsv")

x0<-unique(cancermine_collated[,c("cancer_normalized","gene_normalized")])
x2<-aggregate(x0$cancer_normalized,list(x0$gene_normalized),paste)
rownames(x2)<-x2[,1]
cancermine<-x2

cancermine$Carcinoma<-rep("",length(cancermine[,1]))
cancermine$Sarcoma<-rep("",length(cancermine[,1]))
cancermine$Myeloma<-rep("",length(cancermine[,1]))
cancermine$Leukemia<-rep("",length(cancermine[,1]))
cancermine$Lymphoma<-rep("",length(cancermine[,1]))
cancermine$Mixed <-rep("",length(cancermine[,1]))

cancermine$Carcinoma[grepl("Carcinoma",cancermine[,2] ,ignore.case = T)]="Carcinoma"
cancermine$Carcinoma[grepl("cancer",cancermine[,2]  ,ignore.case = T)]="Carcinoma"
cancermine$Sarcoma[grepl("Sarcoma",cancermine[,2]  ,ignore.case = T)]="Sarcoma"
cancermine$Sarcoma[grepl("astrocytoma",cancermine[,2]  ,ignore.case = T)]="Sarcoma"
cancermine$Sarcoma[grepl("Glioma",cancermine[,2]  ,ignore.case = T)]="Sarcoma"
cancermine$Sarcoma[grepl("Mesenchymous",cancermine[,2] ,ignore.case = T )]="Sarcoma"
cancermine$Myeloma[grepl("Myeloma",cancermine[,2]  ,ignore.case = T)]="Myeloma"

cancermine$Leukemia[grepl("Leukemia",cancermine[,2] ,ignore.case = T )]="Leukemia"
cancermine$Lymphoma[grepl("Lymphoma",cancermine[,2]  ,ignore.case = T)]="Lymphoma"
cancermine$Mixed[grepl("teratocarcinoma",cancermine[,2]  ,ignore.case = T)]="Mixed"
cancermine$Mixed[grepl("carcinosarcoma",cancermine[,2] ,ignore.case = T )]="Mixed"
cancermine$Mixed[grepl("mixed mesodermal tumor",cancermine[,2] ,ignore.case = T )]="Mixed"
cancermine$Mixed[grepl("adenosquamous carcinoma",cancermine[,2] ,ignore.case = T )]="Mixed"

rownames(cancermine)<-cancermine[,1]
cancermine$cancer_types<-rep("other",length(cancermine[,1]))
cancermine$cancer_types <- apply(cancermine[,3:8] , 1 , paste , collapse = "" )
cancermine$cancer_count <- apply(cancermine[,3:8] , 1 ,function(x)length(x[x!=""]) )

cancermine$cancer_types2<-cancermine$cancer_types
cancermine$cancer_types2[cancermine$cancer_count>1]="mixed"
cancermine$cancer_types3<-cancermine$cancer_count

cancermine$cancer_types3[cancermine$cancer_count>0]="less"
cancermine$cancer_types3[cancermine$cancer_count>2]="more"
cancermine$cancer_count2<-cancermine$cancer_count
cancermine$cancer_count2[cancermine$cancer_count>2]=">2"
#cancermine$cancer_count2[is.na(cancermine$cancer_count2)]="0"
genes_infor<-cbind(genes_infor,cancermine[genes_infor$hgnc_symbol,])
genes_infor$cancer_count2<-factor(genes_infor$cancer_count2,levels=c("0","1","2",">2"))





col_x=c("hg38_3UTR_length_log")

hist_list<-c("Carcinoma","Sarcoma","Leukemia","Lymphoma","Myeloma")
p<-plot_histoligical(genes_infor,"hg38_3UTR_length_log",comps=pairs_list(hist_list),"Log10(Length of 3'UTR)")
p1_hist<-plot_histoligical(genes_infor,"hg38_3UTR_length_log",comps=pairs_list(hist_list)[1:3],"Log10(Length of 3'UTR)")



overlap <- read.delim("./data/overlap.tsv")# TSG enriched Go terms
genesets <- getGmt("./data/c5.go.mf.v7.2.symbols.gmt",geneIdType=SymbolIdentifier())

kegg=geneIds(genesets)
kegg_data<-ldply (kegg, data.frame)
colnames(kegg_data)<-c("term","hgnc_symbol")
go_names<-colnames(overlap)[grepl("GO_",colnames(overlap))]

kegg_data_select=kegg_data[kegg_data$term %in% go_names,]
kegg_gene_infor=merge(kegg_data_select,genes_infor)
kegg_gene_infor$term<-gsub("GO_","",kegg_gene_infor$term)
temp_kegg<-aggregate(.~term+Gene.Type,data=kegg_gene_infor[,c("term","Gene.Type","hg38_3UTR_length_log")],median)
temp_kegg_tsg=temp_kegg[temp_kegg$Gene.Type=="TSG",]
temp_kegg_non_cancer=temp_kegg[temp_kegg$Gene.Type=="Non-Cancer",]
rownames(temp_kegg_non_cancer)<-temp_kegg_non_cancer[,1]
diff<-temp_kegg_tsg[,3]-temp_kegg_non_cancer[temp_kegg_tsg[,1],3]
names(diff)<-temp_kegg_tsg[,1]
kegg_gene_infor$term<-factor(kegg_gene_infor$term,levels = names(sort(diff)))


P_GO<-plotSubSig(kegg_gene_infor[kegg_gene_infor$Gene.Type!="Oncogene",],"term", "hg38_3UTR_length_log","Gene.Type","", cols1 = c("#1874CD","#EE2C2C"),"TSG enriched Go terms","Log10(Length of 3'UTR)",0.4) + labs(fill = "Gene.Type")
#+coord_flip()


pall0<-ggarrange(p_utr3_count,p_max_3UTR,p_aver_3UTR,nrow=1,ncol=3,labels = letters[1:3], font.label = list(size = 18, color = "black", face = "bold"),legend = "none")
pall00<-ggarrange(p_age_3utr1,p_age_3utr0,nrow=1,ncol=2,labels = letters[6:7],font.label = list(size = 18, color = "black", face = "bold"),legend = "none")
pall1<-ggarrange(pall0,p0_3utr,p1_3utr,nrow=1,ncol=3,widths = c(3,1,1),labels = c("",letters[4:5]), font.label = list(size = 18, color = "black", face = "bold"),legend = "none")
pall2<-ggarrange(pall00,p1_hist$p1,p1_hist$p2, nrow=1,ncol=3,widths = c(2,1,1),labels = c("",letters[8:9]), font.label = list(size = 18, color = "black", face = "bold"),legend = "none")
pall3<-ggarrange(pall1,pall2,nrow = 2,ncol = 1,labels=c("",""), font.label = list(size = 18, color = "black", face = "bold"))
pall4<-ggarrange(pall3,P_GO,nrow = 1,ncol = 2,widths = c(2,1),labels=c("",letters[10]), font.label = list(size = 18, color = "black", face = "bold"),common.legend = T,legend = "bottom")
tiff("./result/figure_1_3UTR.tiff",width = 15,height = 12.5,res=300,units="in",compression = "lzw")
print(pall4)
dev.off()
















