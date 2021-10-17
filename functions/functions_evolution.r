options(stringsAsFactors = FALSE)
library(rstatix)
library(ggplotify)
library("ggExtra")
library(ggsignif)
library("cowplot")
library("gridExtra")
library("RColorBrewer")
library("ggpubr")
library(ggplot2)
library(stringr)
require(plyr)
library(metacor)
library(foreach)
library(doParallel)
library(meta)
library(parallel)
library(GenomicRanges)
library(reshape2)
library(Repitools)
library(EnrichedHeatmap)
library(evobiR)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(proto)
library(GGally)
require(plyr)
library(png)
library(grid)
library(cocor)
library(metacor)
library(meta)
library(GSEABase)
#install.packages("extrafont")
library(extrafont)
#font_import()
### get the expression wildth

get_cv<-function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}
#get_cv_log<-function(x){sqrt(exp(sd(x,na.rm =T)^2)-1)}
get_cv_log<-function(x){
  
  x=2^x-1##get the original TPM value
  return(sd(x,na.rm=T)/mean(x,na.rm=T))
}
get_TSEI_log2<-function(data_x){
  ##transform the log(TPM+1) back to TPM
  data_x=2^data_x-1
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}


get_TSEI_log2_array<-function(data_x){
  ##transform the log(TPM+1) back to TPM
  data_x=2^data_x
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}


get_TSEI_pre<-function(data_x){
  sum(1-2^(data_x-max(data_x)))/(length(data_x)-1)
}


get_TSEI<-function(data_x){
  sum(1-data_x/max(data_x))/(length(data_x)-1)
}

countFunction <- function(x){
  return(data.frame(y=0,label=round(length(x),2)))}

get_substring<-function(x,sp,n){
  x_1=strsplit(x,sp,fixed = T)
  x_2=do.call(rbind.data.frame, x_1)
  return(as.character(x_2[,n]))
}
theme_hd<-function(base_size = 14, base_family = "Helvetica"){
  theme_classic(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    axis.line = element_line(colour = "black",size = 1),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    strip.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}


theme_hd_minimal<-function(base_size = 14, base_family = "Helvetica"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    axis.line = element_line(colour = "black",size = 1),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    strip.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}

theme_hd_minimal2<-function(base_size = 14, base_family = "Helvetica"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%  theme(plot.title = element_text(lineheight=.8, size=14, face="bold"),
                                                                                    legend.background = element_rect(colour = "white"),
                                                                                    legend.key = element_rect(colour = "white"),legend.position = "bottom",
                                                                                    axis.text = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    axis.title = element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.text =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    legend.title =  element_text( size=14,face="bold",vjust = 0.5),
                                                                                    plot.margin = unit(c(1,1,1,1), "cm"))}



get_main_transcript<-function(get_gene_pos_temp){
  # get the transcript which is principal1
  main_trans<-subset( get_gene_pos_temp,APPRIS.annotation=="principal1")
  genes=unique(get_gene_pos_temp$ensembl_gene_id)
  genes=genes[! genes %in% main_trans$ensembl_gene_id] ## genes which have no principal transcript
  get_gene_pos_temp=get_gene_pos_temp[!get_gene_pos_temp$ensembl_gene_id  %in% main_trans$ensembl_gene_id,]
  get_gene_pos_temp<-rbind(main_trans,get_gene_pos_temp) # genes with  principal1 and genes that do not have principal1
  max_trans_length<-aggregate(get_gene_pos_temp$trans_length,list(get_gene_pos_temp$ensembl_gene_id),function(x){max(x,na.rm = T)}) # get the max transcript length
  f_select=paste(get_gene_pos_temp$ensembl_gene_id,get_gene_pos_temp$trans_length,sep="_") %in% paste(max_trans_length[,1],max_trans_length[,2],sep="_")
  select_gene_trans<-get_gene_pos_temp[f_select,]
  # for genes with transcript which have the same length randomly select one
  return(select_gene_trans[!duplicated(select_gene_trans$ensembl_gene_id),])
}




plot_age_trans_count<-function(all_data_age_x,types_x){
  sts <- boxplot.stats(all_data_age_x$transcript_count)$stats
  p=ggplot(all_data_age_x, aes(x=factor(Age_class), y=transcript_count, fill=Gene.Type)) + 
    geom_boxplot(colour = "#E5E5E5", outlier.color = "#B3B3B3",outlier.shape = NA) + 
    coord_cartesian(ylim = c(sts*1.3,sts/1.3))+
    scale_fill_manual(values=c("#1874CD","#EE2C2C"), name = "Gene Group") +ggtitle(types_x)+
    theme_bw()+
    theme(plot.title = element_text(lineheight=.8, size=14, face="bold.italic"),
          axis.line.x = element_line(colour = "black",size = 0.7,lineend = 2),
          axis.line.y = element_line(colour = "black",size = 0.7),
          axis.title.x = element_text( size=17, face="bold"),
          axis.title.y = element_text( size=17, face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x=element_text( size=12,face="bold",vjust = 0.3),
          axis.text.y=element_text( size=12,face="bold",vjust = 0.3,angle = 90),
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(colour = "white"),
          plot.margin = unit(c(0.5,1,0.5,0.3), "cm")) + labs(x = "Origin Time", y = "transcript count")
  print(p)
}



#####################33 Only calculated one time if two TSG or non-cancer genes are paralogs    ; calculated the mean transcript count for the paralogs 
calculate_paralog<-function(data_x,all_genes,trans_count){
  paralog_transcount<-c()
  genes=unique(data_x$ensembl_gene_id)
  gene_select=c()
  gene_dup=c()
  paralog_unique=c()
  for(gene in genes){
    if(gene %in% gene_dup){next}
    para_gene=data_x[data_x[,1]==gene,]
    types_gene=all_genes[gene,"Gene.Type"]
    temp_cancer=c(gene,para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type==types_gene]])
    temp_non_cancer=para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type=="Non-Cancer"]]
    temp_other_cancer=para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[!all_genes$Gene.Type %in% c("Non-Cancer",types_gene)]]
    temp_cancer_count=c(mean(trans_count$all_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$pro_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$non_pro_trans[temp_cancer,"transcript_count"],rm.na=T),
                        mean(trans_count$non_pro_trans[temp_cancer,"Age_class"],rm.na=T))
    temp_non_cancer_count=c(mean(trans_count$all_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$pro_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$non_pro_trans[temp_non_cancer,"transcript_count"],rm.na=T),
                            mean(trans_count$non_pro_trans[temp_non_cancer,"Age_class"],rm.na=T))
    temp_other_cancer_count=c(mean(trans_count$all_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$pro_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$non_pro_trans[temp_other_cancer,"transcript_count"],rm.na=T),
                              mean(trans_count$non_pro_trans[temp_other_cancer,"Age_class"],rm.na=T))
    paralog_transcount<-rbind(paralog_transcount,c(length(temp_cancer),temp_cancer_count,
                                                   length(temp_non_cancer),temp_non_cancer_count,
                                                   length(temp_other_cancer),temp_other_cancer_count))
    gene_dup=c(gene_dup,temp_cancer)
    gene_select=c(gene_select,gene)
    paralog_unique<-rbind(paralog_unique,para_gene)
  }
  rownames(paralog_transcount)<-gene_select
  colnames(paralog_transcount)<-paste(rep(c("cancer","non_cancer_paralog","other_patalog"),each=5),c("length","all","protein","non_protein","Age_class"),sep="_")
  return(list(a=as.data.frame(paralog_transcount),b=paralog_unique))
}

calculate_paralog_2<-function(data_x,all_genes){
  paralog_transcount<-c()
  genes=unique(data_x$ensembl_gene_id)
  gene_select=c()
  gene_dup=c()
  paralog_unique=c()
  for(gene in genes){
    if(gene %in% gene_dup){next}
    para_gene=data_x[data_x[,1]==gene,]
    types_gene=all_genes[gene,"Gene.Type"]
    temp_cancer=c(gene,para_gene[,2][para_gene[,2] %in% all_genes$Gene.ID[all_genes$Gene.Type==types_gene]])
    gene_dup=c(gene_dup,temp_cancer)
    paralog_unique<-rbind(paralog_unique,para_gene)
  }
  return(paralog_unique)
}


get_paralog_same_age_type<-function(x_paralog,all_gene_trans_infor){
  genes<-unique(x_paralog[,1])
  paralog_pairs<-c()
  for(gene in genes){
    
    temp_p<-x_paralog[x_paralog[,1]==gene,]
    temp_infor<-all_gene_trans_infor[unique(c(gene,temp_p[,2])),]
    temp_age=temp_infor[which(temp_infor$age_type==temp_infor[gene,"age_type"]),]
    if(length(temp_age[,1])==0){next}
    trans_count=aggregate(temp_age$all_trans_count,list(temp_age$Gene.Type),median)
    gene_count=aggregate(temp_age$all_trans_count,list(temp_age$Gene.Type),length)
    pro_count=aggregate(temp_age$pro_count,list(temp_age$Gene.Type),median)
    non_pro_count=aggregate(temp_age$non_pro_trans_count,list(temp_age$Gene.Type),median)
    promoter_cpg_oe=aggregate(temp_age$promoter_cpg_OE,list(temp_age$Gene.Type),median)
    temp_result=cbind(gene,pro_count[,1],temp_infor[gene,"age_type"],gene_count[,2],trans_count[,2], pro_count[,2],non_pro_count[,2],promoter_cpg_oe[,2])
    paralog_pairs<-rbind(paralog_pairs,temp_result)
  }
  rownames(paralog_pairs)<-paralog_pairs[,1]
  colnames(paralog_pairs)<-c("gene","Gene.Type","age_type","Gene.Type.count","all_trans_count","pro_trans_count","non_pro_trans_count","promoter_cpg_OE")
  return(as.data.frame(paralog_pairs))
}





plot_promoter_data<-function(genes_infor,title_x,y_name){
  TSG_non_cancer<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  TSG_non_cancer$HK_types=factor(paste(TSG_non_cancer$HK,TSG_non_cancer$Gene.Type,sep=":"),levels = c("HK_gene:Non-Cancer","non_HK_gene:TSG","HK_gene:TSG","non_HK_gene:Non-Cancer"))
  comp3=list(c("1_13", "13_20"),c("1_13", "20_23"),c("1_13", "23+"),c("13_20", "20_23"),c("13_20", "23+"),c("20_23", "23+"))
  x=quantile(TSG_non_cancer[,y_name], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  s1=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(age_type) ,y=get(y_name),color=HK))  + 
    facet_wrap(~Gene.Type+HK,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                                   textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*3,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*9))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_type")
  
  s2=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(Age_class) ,y=get(y_name),color=I("white"),fill=HK))  + 
    facet_wrap(~Gene.Type+HK,nrow=1) +
    geom_boxplot(alpha=I(0.7),position = position_dodge(width=0.9))+scale_y_continuous(limits =x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_class")+ggtitle("")
  p_all<-ggarrange(s1, s2, nrow = 2, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )
  print(p_all)
}


plot_promoter_data_2<-function(genes_infor,title_x,y_name){
  TSG_non_cancer<-genes_infor[genes_infor$Gene.Type!="Oncogene",]
  comp3=list(c("Non-Cancer", "TSG"))
  x=quantile(TSG_non_cancer[,y_name], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  s1=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=Gene.Type ,y=get(y_name),col=Gene.Type))  + 
    facet_wrap(~age_type,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                               textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*6,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*7))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_type")
  
  s2=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=Gene.Type ,y=get(y_name),col=Gene.Type))  + 
    facet_wrap(~HK+age_type,nrow=1) + geom_signif(comparisons = comp3,map_signif_level=TRUE,
                                                  textsize=3,tip_length = 2,col="black", test='wilcox.test', y_position=seq(x[2]+steps*3,x[2]+steps*9,steps))+
    geom_point(size=I(0.5),position=position_jitterdodge(dodge.width=0.9)) +
    geom_boxplot(alpha=I(0.7),outlier.colour = NA, 
                 position = position_dodge(width=0.9))+scale_y_continuous(limits =x+c(0,steps*7))+ggtitle(title_x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=6))+xlab("Age_type")
  s3=ggplot(TSG_non_cancer[!is.na(TSG_non_cancer$Age_class),],aes(x=as.factor(Age_class) ,y=get(y_name),color=I("white"),fill=Gene.Type))  + 
    facet_wrap(~Gene.Type,nrow=2) +
    geom_boxplot(alpha=I(0.7),position = position_dodge(width=0.9))+scale_y_continuous(limits =x)+ylab(y_name)+theme(axis.text.x = element_text(angle = 90))+xlab("Age_class")+ggtitle("")
    p_all<-ggarrange(s1, s2, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom"  )
    p_all1<-ggarrange(p_all, s3, ncol  = 2, common.legend = TRUE, legend = "bottom" )
  print(p_all1)
}




firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



plot_boxplot<-function( all_data_tsg,y_lab,title_x){
  comp=list(c("low","median"),c("low","high"),c("median","high"))
  x=quantile(all_data_tsg[,y_lab], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  p=ggplot(data = all_data_tsg,aes(y=get(y_lab),x=promoter_cpg_OE_category,fill=promoter_cpg_OE_category,col=I("grey")))+
    geom_boxplot(outlier.colour = NA)+facet_wrap(~cancer,nrow=1)
  s3=p+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=25),axis.title =element_text(size=20))+ 
    scale_y_continuous(limits =x+c(0,steps*5))+
    geom_signif(comparisons = comp,map_signif_level=TRUE,col="black",tip_length = 4, test='wilcox.test', y_position=seq(x[2]+steps*2,x[2]+steps*9,steps))+ggtitle(title_x)
  
  return(s3+ylab(y_lab))
}

plot_boxplot_2<-function( all_data_tsg,y_lab,title_x){
  comp=list(c("low","median"),c("low","high"),c("median","high"))
  x=quantile(all_data_tsg[,y_lab], c(0.1, 0.9),na.rm=T)
  steps=(x[2]-x[1])/5
  p=ggplot(data = all_data_tsg,aes(y=get(y_lab),x=promoter_cpg_OE_category,fill=Gene.Type,col=I("grey")))+
    geom_boxplot(outlier.colour = NA)+facet_wrap(~cancer,nrow=1)
  s3=p+theme(axis.text.x = element_text(angle = 90),strip.text = element_text(size=25),axis.title =element_text(size=20))+ 
    scale_y_continuous(limits =x+c(0,steps*5))+
    geom_signif(comparisons = comp,map_signif_level=TRUE,col="black",tip_length = 4, test='wilcox.test', y_position=seq(x[2]+steps*2,x[2]+steps*9,steps))+ggtitle(title_x)
  
  return(s3+ylab(y_lab))
}






### delete low expression
get_cor_for_each_cancer<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    print(cancer)
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff_n=quantile(temp_data$norm_exp,0.2,na.rm=T)
    quantile_cutoff_t=quantile(temp_data$tumor_exp,0.2,na.rm=T)
    if (cancer %in% c("CRC_GSE137327_RNA","Lung_IMA_GSE86958_RNA")){
      quantile_cutoff_n=1
      quantile_cutoff_t=1   
    }
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff_n & temp_data$tumor_exp< quantile_cutoff_t )),]
     temp_data$diff_exp=abs(temp_data$diff_exp)
        temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }

    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}

get_cor_for_each_cancer_tcga<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    #print(cancer)
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1 
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff & temp_data$tumor_exp< quantile_cutoff )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }
    
    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}

get_filter_low_expression_DW_tcga<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff=1
    print(quantile_cutoff)   
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff & temp_data$tumor_exp< quantile_cutoff )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    result=rbind(result,temp_data)
  }
  return(result[result$Gene.Type=="TSG",])
}


get_filter_low_expression_DW<-function(all_data_tsg){
  cancers=unique(all_data_tsg$cancer)
  result=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    quantile_cutoff_n=quantile(temp_data$norm_exp,0.2,na.rm=T)
    quantile_cutoff_t=quantile(temp_data$tumor_exp,0.2,na.rm=T)
    if (cancer %in% c("CRC_GSE137327_RNA","Lung_IMA_GSE86958_RNA")){
      quantile_cutoff_n=1
      quantile_cutoff_t=1   
    }
    temp_data=temp_data[which( is.finite(temp_data$diff_exp) & temp_data$diff_exp<0 & !(temp_data$norm_exp< quantile_cutoff_n & temp_data$tumor_exp< quantile_cutoff_t )),]
    temp_data$diff_exp=abs(temp_data$diff_exp)
    result=rbind(result,temp_data)
     }
  return(result[result$Gene.Type=="TSG",])
}

get_cor_for_each_cancer0<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    temp_cor_non_cancer=cor.test(temp_data[temp_data$Gene.Type=="Non-Cancer",name1],temp_data[temp_data$Gene.Type=="Non-Cancer",name2],method = method_x)
    temp_cor_TSG=cor.test(temp_data[temp_data$Gene.Type=="TSG",name1],temp_data[temp_data$Gene.Type=="TSG",name2],method = method_x)
    len_non_cancer=length(temp_data[temp_data$Gene.Type=="Non-Cancer",name1])
    len_TSG=length(temp_data[temp_data$Gene.Type=="TSG",name1])
    if(method_x=="kendall"){
      r1=sin(pi*0.5*temp_cor_non_cancer$estimate)
      r2=sin(pi*0.5*temp_cor_TSG$estimate)      
    }else{
      r1=temp_cor_non_cancer$estimate
      r2=temp_cor_TSG$estimate    
    }

    if(r1>r2){temp="greater"
    temp_x="<"}else{temp="less"
    temp_x=">"}
    cocor_x=cocor.indep.groups(r1.jk=r1, r2.hm=r2, n1=len_non_cancer, n2=len_TSG, alternative=temp, alpha=0.05, conf.level=0.95, null.value=0)
    result_cor<-rbind(result_cor,c(cancer,"Non_Cancer",round(temp_cor_non_cancer$estimate,3),temp_cor_non_cancer$p.value,cocor_x@fisher1925$p.value))
    result_cor<-rbind(result_cor,c(cancer,"TSG",round(temp_cor_TSG$estimate,3),temp_cor_TSG$p.value,cocor_x@fisher1925$p.value))
  }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue","Diff_pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[result_cor$Pvalue==0]=2.2e-16
  result_cor$Pvalue=-log10(as.numeric(result_cor$Pvalue))
  result_cor$Diff_pvalue=as.numeric(result_cor$Diff_pvalue)
  return(as.data.frame(result_cor))
}


get_cor_for_each_cancer_all<-function(all_data_tsg,name1,name2,method_x){
  cancers=unique(all_data_tsg$cancer)
  result_cor=c()
  for(cancer in cancers){
    temp_data=all_data_tsg[all_data_tsg$cancer==cancer,]
    temp_cor_all=cor.test(temp_data[,name1],temp_data[,name2],method = method_x)
    result_cor<-rbind(result_cor,c(cancer,"All_genes",round(temp_cor_all$estimate,3),formatC(temp_cor_all$p.value, format = "e", digits = 2)))
 }
  colnames(result_cor)<-c("Cancer","Gene.Type","cor","Pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  
  if(result_cor$cor[1]<0){order_cancer=result_cor[order(result_cor$cor,decreasing = T),1]}else{order_cancer=result_cor[order(result_cor$cor,decreasing = F),1]}
  
  result_cor$Cancer=factor(result_cor$Cancer,levels=order_cancer)
  return(as.data.frame(result_cor))
}

get_cor_for_each_age<-function(all_data_tsg,name1,name2,ylab_x,method_x){
  cancers=unique(all_data_tsg$cancer)
  ages=unique(all_data_tsg[!is.na(all_data_tsg[,ylab_x]),ylab_x])
  result_cor=c()
  for(cancer in cancers){
    for(age in ages){
      temp_data=all_data_tsg[all_data_tsg[,ylab_x]==age &all_data_tsg$cancer==cancer ,]
      temp_cor_all=cor.test(temp_data[,name1],temp_data[,name2],method = method_x)
      result_cor<-rbind(result_cor,c(cancer,age,"All_genes",round(temp_cor_all$estimate,3),formatC(temp_cor_all$p.value, format = "e", digits = 2)))
      
    }
  }
  colnames(result_cor)<-c("Cancer","age","Gene.Type","cor","Pvalue")
  result_cor=as.data.frame(result_cor)
  result_cor$cor=as.numeric(result_cor$cor)
  result_cor$Pvalue[as.numeric(result_cor$Pvalue)<2.2e-16]=2.2e-16
  
  return(as.data.frame(result_cor))
}



plot_cor<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  #print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                 ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                             yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))
  
  cor_cpg_exp$Significance<- ifelse(round(10^-cor_cpg_exp$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  
   cor_plot1=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue,shape=Significance))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))+
     scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  
   #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot1=cor_plot1,cor_plot_matrix=cor_cpg_exp))
}



plot_cor_all<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer_all(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="All_genes",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))
  
  cor_cpg_exp$Significance<- ifelse(round(10^-cor_cpg_exp$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  
  cor_plot1=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                     yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue,shape=Significance))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(labels=c("Non-Cancer","TSG"), values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")+xlab(firstup(method_x))+
    scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot1=cor_plot1,cor_plot_matrix=cor_cpg_exp))
}


plot_cor_1<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("-log10(0.05)","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1.3,5,10,15),guide="legend")+xlab(firstup(method_x))
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot_matrix=cor_cpg_exp))
}


plot_cor_no_order<-function(all_data_tsg,name1,name2,title_x,method_x){
  cor_cpg_exp<-get_cor_for_each_cancer0(all_data_tsg,name1,name2,method_x)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>15]=15.7
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  order_cancer=levels(all_data_tsg$cancer)
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  print(cor_cpg_exp)
  cor_cpg_exp_tsg=cor_cpg_exp[cor_cpg_exp$Gene.Type=="TSG",]
  Sig <- ifelse(cor_cpg_exp_tsg$Diff_pvalue>0.05," ",ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.01, "*" , 
                                                            ifelse(cor_cpg_exp_tsg$Diff_pvalue >0.001, "**" , "***")))
  
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=cor))+geom_segment(data=cor_matrix,aes(x=TSG,        xend=Non_Cancer, y=factor(Cancer), 
                                                                                    yend=factor(Cancer)),colour="#BFEFFF", size = 1)+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+
    geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label=Sig)+ geom_point(data=cor_cpg_exp,aes(y=Cancer,x=cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("-log10(0.05)","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1.3,5,10,15),guide="legend")+xlab(firstup(method_x))
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(list(cor_plot=cor_plot,cor_plot_matrix=cor_cpg_exp))
}

plot_cor_TSG<-function(all_data_tsg,name1,name2,title_x){
  cor_cpg_exp<-get_cor_for_each_cancer(all_data_tsg,name1,name2)
  cor_cpg_exp$Pvalue[cor_cpg_exp$Pvalue>30]=29
  cor_matrix=reshape2::dcast(cor_cpg_exp[,1:3],Cancer~Gene.Type)
  cor_matrix$mean_cor=(cor_matrix$Non_Cancer+cor_matrix$TSG)/2
  if(cor_matrix$mean_cor[1]>0){order_cancer=cor_matrix[order(cor_matrix$Non_Cancer),1]}else{order_cancer=cor_matrix[order(cor_matrix$TSG),1]}
  
  cor_cpg_exp$Cancer=factor(cor_cpg_exp$Cancer,levels=order_cancer)
  cor_matrix$Cancer=factor(cor_matrix$Cancer,levels=order_cancer)
  cor_plot=ggplot(cor_cpg_exp,aes(y=Cancer,x=Pearson_Cor))+geom_segment(data=cor_matrix,aes(x=TSG,xend=Non_Cancer, y=factor(Cancer), yend=factor(Cancer)),colour="#BFEFFF", size = 1)+
    scale_fill_manual(values=c("#1874CD","#EE2C2C"))+ geom_text(data=cor_matrix,aes(x=mean_cor,y=factor(Cancer)),label="***")+
    geom_point(data=cor_cpg_exp,aes(y=Cancer,x=Pearson_Cor,col=Gene.Type,size=Pvalue))+labs(size = "-Log10(P-value)")+ggtitle(title_x)+
    scale_color_manual(values=c("#1874CD","#EE2C2C"))+scale_size_area(labels=c("1","5","10",">15"),max_size=7,limits=c(0,16),breaks=c(1,5,10,15),guide="legend")
  #theme(legend.position = "bottom",legend.box = "vertical")+xlab("Pearson_cor\nexpression in normal sample VS Promoter CpG OE")
  return(cor_plot)
}




plot_non_paired_paralog<-function(TSG_paralog_non_cancer,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non_Cancer_paralogs"))
  p1=ggplot(data = TSG_paralog_non_cancer[!is.na(TSG_paralog_non_cancer[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    ggtitle("Paralogs")+ theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}


plot_TSG_non_cancer<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=-0.5,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  return(p1+ylab(y_lab))
}

plot_TSG_non_cancer_withlegend<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=0,label=round(length(x),2)))}
  
  comp=list(c("TSG", "Non-Cancer"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = TRUE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=function(p)sprintf("p = %.2g", p),col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"),axis.text.x =element_blank(),axis.title.x = element_blank())
  p1=p1+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}

plot_all_genes<-function(genes_infor,y_name,y_lab){
  countFunction <- function(x){
    return(data.frame(y=0,label=round(length(x),2)))}
  genes_infor$Gene.Type=factor(genes_infor$Gene.Type,levels = c("Non-Cancer","TSG","Oncogene"))
  comp=list(c("TSG", "Non-Cancer"),c("TSG","Oncogene"),c("Non-Cancer","Oncogene"))
  p1=ggplot(data = genes_infor[!is.na(genes_infor[,y_name]),],aes(x=Gene.Type,y=get(y_name),fill=Gene.Type))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=function(p)sprintf("p = %.2g", p),col="black",step_increase=0.1,textsize = 4.8,fontface="bold") +scale_fill_manual(values=c("#1874CD","#2E8B57", "#EE2C2C"))+theme_hd()+
    theme(plot.title = element_text(face = "bold"))
  p1=p1+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")+ylab(y_name)
  return(p1+ylab(y_lab))
}


plotPairedParalog<-function(TSG_paralog_age,col.x,lab.x){
  
  tsg_genes=TSG_paralog_age[TSG_paralog_age[,2]=="TSG",]
  non_cancer_genes=TSG_paralog_age[TSG_paralog_age[,2]=="Non-Cancer",]
  overlap_gene=as.character(tsg_genes[as.character(tsg_genes[,1]) %in% as.character(non_cancer_genes[,1]),1])
  
  data_x=data.frame(TSG=as.numeric(as.character(tsg_genes[overlap_gene,col.x])),Non_Cancer_Paralogs=as.numeric(as.character(non_cancer_genes[overlap_gene,col.x])),age=as.numeric(as.character(non_cancer_genes[overlap_gene,]$Age_class)))
  
  p2=ggpaired(data_x[!(is.na(data_x$TSG) | is.na(data_x$Non_Cancer_Paralogs)),], cond2 = "TSG", cond1 = "Non_Cancer_Paralogs",
              color  = "condition", line.color = "gray78", line.size = 0.1,palette = "npg",size=1,width = 0.5)+ylab(lab.x)+theme_hd()+ theme(plot.title = element_text(hjust = 0.5,face = "plain"),axis.text.x  = element_blank(),axis.title.x = element_blank())+
    stat_compare_means(method="wilcox.test", paired=TRUE, aes(label = paste0("p = ", ..p.format..)),size=4.8,fontface="bold")+ggtitle("Paralogs with the same age")+scale_color_manual("Condition",values=c("#1874CD","#EE2C2C"))
  p2=p2+theme(axis.title.y = element_text(size = 14),plot.title = element_text(size=14,face="bold"))
  return(p2)
}

pairs_list<-function(x){
  re_all<-list()
  
  for(i in 1:(length(x)-1)){
    for(j in (i+1):length(x)){
    re_all[[paste(i,j)]]<-c(x[i],x[j])      
    }
  }
  return(re_all)
}


plotAgelength<-function(genes_infor,y.name,y_lab,pos,y_lab_2){
  genes_infor=genes_infor[!is.na(genes_infor$Age_class),]
  genes_infor$age_type=genes_infor[,y_lab_2]
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.5
  label_x<-unique(genes_infor$age_type)
  comp3=pairs_list(label_x)
  s3=ggplot(genes_infor[!is.na(genes_infor$Age_class),],aes(x=as.factor(age_type) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(age_type=label_x,y.n = rep(in.x1,length(label_x)))
  colnames(label.df)<-c("age_type",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(age_type = x, x= y)
  colnames(arc.df)<-c("age_type",y.name)
  t.p<-c()
  for(i in label.df$age_type){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$age_type==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$age_type==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(age_type+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+4, get(y.name)+in.x))+xlab("Gene Age (million years)")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()+theme(axis.text.x =   element_text(angle = 30,hjust = 1))
  return(s3)
         #+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}



plotAgelength5utrratio<-function(genes_infor,y.name,y_lab,pos){
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.996,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.996,na.rm = T)*1.5
  comp3=list(c("1_13", "13_20"),c("1_13", "20_23"),c("1_13", "23+"),c("13_20", "20_23"),c("13_20", "23+"),c("20_23", "23+"))
  s3=ggplot(genes_infor[!is.na(genes_infor$Age_class),],aes(x=as.factor(age_type) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(age_type= c("1_13", "13_20","20_23", "23+"),y.n = rep(in.x1,4))
  colnames(label.df)<-c("age_type",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(age_type = x, x= y)
  colnames(arc.df)<-c("age_type",y.name)
  t.p<-c()
  for(i in label.df$age_type){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$age_type==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$age_type==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(age_type+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(age_type+4, get(y.name)+in.x))+xlab("Gene Age class")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()
  return(s3+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}


plotForest<-function(all_data,x.name,y.name,title.x,squaresize,spacing,colgap,method_x="pearson"){
  cor_plot1<-get_cor_for_each_cancer_all(all_data,x.name,y.name,method_x)
  x<-table(all_data$cancer)
  cols.x<-rep("black",length(cor_plot1$Cancer))
  cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
  m1<<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
  p.ff<-grid.grabExpr(print(forest(m1,fontsize = 11,squaresize=squaresize,plotwidth="5cm",colgap = colgap,spacing = spacing,layout ="JAMA" ,ff.fixed="bold",fs.fixed = 13,col.square=cols.x,xlab=title.x,ff.study="bold",fs.heading  = 13,ff.test.overall="bold",fs.xlab = 12,ff.xlab = "bold")))
  print(forest(m1,fontsize = 11,squaresize=0.9,plotwidth="5cm",colgap = colgap,spacing = spacing,layout ="JAMA" ,ff.fixed="bold",fs.fixed = 13,col.square=cols.x,xlab=title.x,ff.study="bold",fs.heading  = 13,ff.test.overall="bold",col.diamond="red",fs.xlab = 12,ff.xlab = "bold"))
  return(p.ff)
}



get_cor_age<-function(all_data,x.name,y.name){
  age<-unique(all_data[!is.na(all_data$Age_class),]$Age_class)
  all_cor<-c()
  for( age_x in age){
      all_data_sub=all_data[which(all_data$Age_class==age_x),]
      cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,"pearson")
      x<-table(all_data_sub$cancer)
      cols.x<-rep("black",length(cor_plot1$Cancer))
      cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
      m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
      all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor$age=factor(all_cor$age,levels = 1:26)
  return(all_cor)
}

get_cor_agetypeprevious<-function(all_data,x.name,y.name,y_lab,method="pearson"){
  all_data$age_type=all_data[,y_lab]
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,"pearson")
    x<-table(all_data_sub$cancer)
    cols.x<-rep("black",length(cor_plot1$Cancer))
    cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
    m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
    print(cor_plot1)
    all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor,p_value=cor_plot1$Pvalue))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
  all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
  all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
  all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))
  
  all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_cor_agetype<-function(all_data,x.name,y.name,y_lab,method="pearson"){
  all_data$age_type=all_data[,y_lab]
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_plot1<-get_cor_for_each_cancer_all(all_data_sub,x.name,y.name,method)
    x<-table(all_data_sub$cancer)
    cols.x<-rep("black",length(cor_plot1$Cancer))
    cols.x[as.numeric(cor_plot1$Pvalue)<0.05]="red"
    m1<-metacor(cor_plot1$cor,as.numeric(x[cor_plot1$Cancer]),studlab = cor_plot1$Cancer)
    print(cor_plot1)
    all_cor<-rbind(all_cor,cbind(cancer=m1$studlab,age_number=m1$n, age=age_x,fix=m1$TE.fixed,random=m1$TE.random,cor= m1$cor,p_value=cor_plot1$Pvalue))
  }
  all_cor<-as.data.frame(all_cor)
  all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
  all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
  all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
  all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))
  
    all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_cor_agetype_single<-function(all_data,x.name,y.name,y_lab){
  all_data$age_type=all_data[,y_lab]
  age<-unique(all_data[!is.na(all_data$age_type),]$age_type)
  all_cor<-c()
  for( age_x in age){
    all_data_sub=all_data[which(all_data$age_type==age_x),]
    cor_x<-cor.test(all_data_sub[,x.name],all_data_sub[,y.name])
    all_cor<-rbind(all_cor,cbind(age=age_x,cor= cor_x$estimate,p_value=cor_x$p.value))
  }

all_cor<-as.data.frame(all_cor)
all_cor[,"cor"]<-as.numeric(all_cor[,"cor"])
all_cor[,"p_value"]<-as.numeric(all_cor[,"p_value"])
all_cor$Significance <- ifelse(all_cor[,"p_value"]>0.05,"NS",ifelse(all_cor[,"p_value"] >0.01, "*" , ifelse(all_cor[,"p_value"] >0.001, "**" , "***")))
all_cor$Significance <- factor(all_cor$Significance ,levels = c("NS","*","**","***"))

all_cor$age=factor(all_cor$age,levels = levels(all_data$age_type))
  return(all_cor)
}


get_diff<-function(data_x,x_lab,y_lab,pos_x){
  re.result<-c()
  cancers=unique(data_x$cancer)
  data_x=data_x[data_x$Gene.Type!="Oncogene",]
  for(cancer in cancers){
    temp=data_x[data_x$cancer==cancer,c(y_lab,x_lab)]
    colnames(temp)<-c("yvalue","x_lab")
    x.p<-wilcox.test(temp$yvalue ~ temp$x_lab)
    temp$x_lab=ifelse(temp$x_lab==pos_x,"pos","neg")
    sum_temp<-temp %>% dplyr::group_by(x_lab) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}


getsigm6a_diff<-function(data_x,m6a_utr,n,y_lab){
  re.result<-c()
  x=table(m6a_utr[m6a_utr$m6a=="+",]$ensembl_gene_id)
  data_x$m6a=0
  data_x$m6a=x[data_x$ensembl_gene_id]
  data_x$m6a[is.na(data_x$m6a)]=0
  data_x$m6a_in<-ifelse(data_x$m6a>=n,"+","no")
  data_x$m6a_in[data_x$m6a==0]<-"-"
  data_x=data_x[which(data_x$m6a_in!="no"),]
  cancers=unique(data_x$cancer)
  for(cancer in cancers){
    temp=data_x[data_x$cancer==cancer,c(y_lab,"m6a_in")]
    colnames(temp)<-c("yvalue","m6a_in")
    temp$m6a_in=ifelse(temp$m6a_in=="+","pos","neg")
    x.p<-wilcox.test(temp$yvalue ~ temp$m6a_in)
    
    sum_temp<-temp %>% dplyr::group_by(m6a_in) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}


getsigm6a_difflen<-function(data_x,n,y_lab,ylab2){
  re.result<-c()
  cancers=unique(data_x$cancer)
  for(cancer in cancers){
    temp=data_x[data_x$cancer==cancer,c(y_lab,ylab2)]
    colnames(temp)<-c("yvalue","m6a_in")
    temp$m6a_in=ifelse(temp$m6a_in=="+","pos","neg")
    x.p<-wilcox.test(temp$yvalue ~ temp$m6a_in)
    print(table(temp$m6a_in))
    sum_temp<-temp %>% dplyr::group_by(m6a_in) %>% dplyr::summarise(dplyr::across("yvalue",list(mean=~mean(.x, na.rm = TRUE),length=length,sd=~sd(.x, na.rm = TRUE))))
    temp.sum<-t(data.frame(as.numeric(c(sum_temp[1,-1],sum_temp[2,-1]))))
    colnames(temp.sum)<-paste(rep(as.character(unlist(sum_temp[,1])),each=3),colnames(sum_temp)[-1],sep = "_")
    re.result<-rbind(re.result,data.frame(cancer=cancer,wilcoxP=x.p$p.value,temp.sum))
  }
  return(re.result)
}




#Potential 3UTR m6A+
plotForesttwogroup<-function(sif_m6a,title_x, colgap,pos,control,spacingx =0.9){
  sif_m6a[,"pos_yvalue_mean"]=round(sif_m6a[,"pos_yvalue_mean"],2)
  sif_m6a[,"neg_yvalue_mean"]=round(sif_m6a[,"neg_yvalue_mean"],2)
  m1 <<- metacont(pos_yvalue_length,  pos_yvalue_mean,pos_yvalue_sd ,neg_yvalue_length,  neg_yvalue_mean,neg_yvalue_sd,studlab=cancer,data =sif_m6a, sm = "SMD",label.e=pos,label.c=control)
  sig1<-sif_m6a[,"wilcoxP"]
  cols.x<-rep("black",length(sif_m6a[,"wilcoxP"]))
  cols.x[as.numeric(sif_m6a[,"wilcoxP"])<0.05]="red"
  cols.x<<-cols.x ### forest only use global varaint
  title_x<<-title_x
  colgap <<- colgap
  spacingx <<-spacingx
  SMD=m1$TE
  names(SMD)<-m1$data$cancer
  p1_3utr_diff <- as.grob(~forest(m1,xlab=title_x,fontsize = 11,squaresize=0.9,spacing = spacingx, col.square =cols.x,colgap = colgap,fs.xlab = 12,ff.xlab = "bold"))
  return(list(p1=p1_3utr_diff,SMD=SMD))
}


getSMDforestcutoff<-function(data_x,m6a_3utr,y_lab,pos,control){
  num<-min(16,length(unique(m6a_3utr[,"cancer"])))
  re_all<-c()
  for(n in 1:num){
    temp_sig<-getsigm6a_diff(data_x,m6a_3utr,n,y_lab)
    print(n)
   m1 <- metacont(pos_yvalue_length,  pos_yvalue_mean,pos_yvalue_sd ,neg_yvalue_length,  neg_yvalue_mean,neg_yvalue_sd,studlab=cancer,data =temp_sig, sm = "SMD",label.e=pos,label.c=control)
   re_all<-rbind(re_all,cbind(cutoff=paste(">=",n),smd=m1$TE,cancer=m1$studlab))
     }

  re_all<-as.data.frame(re_all)
  re_all[,"cutoff"]<-factor(re_all[,"cutoff"],levels =paste(">=",1:num) )
  re_all[,"smd"]<-as.numeric(re_all[,"smd"])
  return(as.data.frame(re_all))
}




plotCorMultiple<-function(all_data_tsg_all_genes,y.name,x1,x1.name,title.x,cols){
  all_norm_plot<-c()
  for(i in 1:length(x1)){
    cor_norm_exp<-get_cor_for_each_cancer_all(all_data_tsg_all_genes,y.name,x1[i],"pearson")
    all_norm_plot<-rbind(all_norm_plot,cbind(type_x=x1.name[i],cor_norm_exp))
  }
  #col=c("#FFC125","#00B2EE","#EE6363", "#458B00", "purple")
  all_norm_plot$Pvalue<-as.numeric(all_norm_plot$Pvalue)
  cor_matrix1= aggregate(.~type_x,all_norm_plot[,-c(2:3)],median)
  print(cor_matrix1)
  cor_matrix2= aggregate(.~type_x,all_norm_plot[,-c(2:3)],mad)
  print(cor_matrix2)
  all_norm_plot$type_x<-factor(all_norm_plot$type_x,levels = x1.name)
  all_norm_plot$Pvalue<--log10(as.numeric(all_norm_plot$Pvalue))
  all_norm_plot$Significance<- ifelse(round(10^-all_norm_plot$Pvalue,5)>0.05,"NS", "Pvalue<0.05")
  if(length(unique(all_norm_plot$Significance))>1){
    cor_plot2= ggplot(all_norm_plot,aes(y=Cancer,x=cor,col=type_x))+  geom_point(aes(size=Pvalue,shape=Significance))+scale_shape_manual("Significance",values=c(1,19),label=c("No","yes"))
  }else{
    cor_plot2= ggplot(all_norm_plot,aes(y=Cancer,x=cor,col=type_x))+  geom_point(aes(size=Pvalue))
    
  }
  cor_plot2=cor_plot2+labs(size = "-Log10(P-value)")+scale_color_manual("",values=cols,label=x1.name)+ggtitle(title.x)+  scale_size_area(labels=c("1","2",">5"),max_size=7,limits=c(0,16),breaks=c(1,2,5))+guides(color=guide_legend(nrow=2,byrow=TRUE),size=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))+xlab("Pearson's r")+geom_hline(yintercept = 0.4,size=1.5)+geom_vline(xintercept = 0,size=1)+theme_hd_minimal2()+ylab("Dataset")+guides(size = guide_legend(order = 1),col = guide_legend(order = 3))
  return(cor_plot2)
}


plotcancercount<-function(genes_infor,y.name,y_lab,pos){
  countFunction <- function(x){
    return(data.frame(y=pos,label=round(length(x),2)))}
  in.x1<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.7
  in.x<-quantile(genes_infor[,y.name],0.8,na.rm = T)*1.5
  comp3=list(c("0", "1"),c("1", "2"),c("2", ">2"))
  s3=ggplot(genes_infor[!is.na(genes_infor$cancer_count2) & genes_infor$Gene.Type!="Oncogene" ,],aes(x=as.factor(cancer_count2) ,y=get(y.name)))  + geom_signif(comparisons = comp3,map_signif_level=TRUE,col="black", test='wilcox.test', step_increase=0.1,tip_length = 0.01) +geom_boxplot(alpha=I(0.7),outlier.colour = NA,aes(fill=Gene.Type),width=0.4)+ylab(y_lab)
  label.df <- data.frame(cancer_count2= c("0", "1","2", ">2"),y.n = rep(in.x1,4))
  colnames(label.df)<-c("cancer_count2",y.name)
  r <- .1
  t <- seq(0, 180, by = 1) * pi / 180
  x <- r * cos(t)
  y <- r*2*in.x1/3 * sin(t)
  arc.df <- data.frame(cancer_count2 = x, x= y)
  colnames(arc.df)<-c("cancer_count2",y.name)
  t.p<-c()
  for(i in label.df$cancer_count2){
    t1<-wilcox.test(genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2==i,y.name],genes_infor[genes_infor$Gene.Type=="Non-Cancer" & genes_infor$cancer_count2==i,y.name])
    t.p<-c(t.p,t1$p.value)
  }
  
  Sig <- ifelse(t.p>0.05,"NS",ifelse(t.p >0.01, "*" , ifelse(t.p >0.001, "**" , "***")))
  
  s3=s3 + geom_text(data = label.df, label = Sig,col=I("blue"))+ geom_line(data = arc.df, aes(cancer_count2+1, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+2, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+3, get(y.name)+in.x))+ geom_line(data = arc.df, aes(cancer_count2+4, get(y.name)+in.x))+xlab("Gene Age class")+scale_fill_manual(values=c("#1874CD","#EE2C2C"))+theme_hd()
  return(s3+stat_summary(fun.data =  countFunction, col="black",geom="text", size = 4.8,fontface="bold"))
}



plotSubSig<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(size_x==0.5){
     stat.test=stat.test %>% filter(p.adj.signif!="ns")
     #print(stat.test)
  }
  if(dim(stat.test)[1]>0){
     p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin,size=size_x, xmax=stat.test$xmax, annotations=stat.test$p.adj.signif, y_position=max(stat.test$y.position),fontface="bold") 
 
  }
  return(p1)
}


plotSubSigOneSide<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black",side){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group,alternative=side)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(size_x==0.5){
    stat.test=stat.test %>% filter(p.adj.signif!="ns")
    #print(stat.test)
  }
  if(dim(stat.test)[1]>0){
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin,size=size_x, xmax=stat.test$xmax, annotations=stat.test$p.adj.signif, y_position=max(stat.test$y.position),step_increase = 0.1,fontface="bold") 
    
  }
  return(p1)
}



plotSubSigall4<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  p<-list()
  cancer=unique(x.data$cancer)
  i=1
  max_x<-max(x.data[,ylab1])+1
  if(length(cancer) %% 4!=0){
    nr<-length(cancer) %/% 4 
  }else{nr<-length(cancer) %/% 4 -1 }
  j=1
  while(i <length(cancer)){
    print(length(cancer))
    if(i%/% 4 !=nr){
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
     p_1<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4,widths = c(1.2,1,1,1))
      
    }else{
      
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+xlab("")
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p_1<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4,widths = c(1.2,1,1,1),common.legend = T,legend = "bottom")
      
    }
    print(i)
   print(p_1)
    p[[j]]=p_1
    i=i+4
    j=j+1
    print(i)
  }
  p_all<-ggarrange(p[[1]],p[[2]],p[[3]],nrow=3,ncol=1,common.legend = T,legend="bottom",heights = c(1,1,1.8))
  return(p_all)  
}

plotSubSigall5<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,size_x=0.5,cols_x="black"){
  p<-list()
  cancer=unique(x.data$cancer)
  i=1
  max_x<-max(x.data[,ylab1])
  if(length(cancer) %% 5!=0){
    nr<-length(cancer) %/% 5 
  }else{nr<-length(cancer) %/% 5 -1 }
  j=1
  while(i <length(cancer)){
    print(length(cancer))
    if(i%/% 5 !=nr){
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p5<-plotSubSig(x.data[x.data$cancer==cancer[i+4],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+4])
      p5=p5+theme(axis.text.x = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),legend.position = "none")+ scale_x_discrete(expand = c(0, 0))+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))
      p_1<-ggarrange(p1,p2,p3,p4,p5,nrow=1,ncol=5,widths = c(1.2,1,1,1,1))
      
    }else{
      
      p1<-plotSubSig(x.data[x.data$cancer==cancer[i],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i])
      p1=p1+xlab("")
      p2<-plotSubSig(x.data[x.data$cancer==cancer[i+1],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+1])
      p2=p2+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p3<-plotSubSig(x.data[x.data$cancer==cancer[i+2],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+2])
      p3=p3+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))
      p4<-plotSubSig(x.data[x.data$cancer==cancer[i+3],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+3])
      p4=p4+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous( limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p5<-plotSubSig(x.data[x.data$cancer==cancer[i+4],],xlab1,ylab1,sub1,title.x,cols1,xlab.name,ylab.name,size_x,cols_x)+ggtitle(cancer[i+4])
      p5=p5+theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")+ scale_y_continuous(limits=c(0,max_x+1),expand = c(0, 0))+ scale_x_discrete(expand = c(0, 0))+xlab("")
      p_1<-ggarrange(p1,p2,p3,p4,p5,nrow=1,ncol=5,widths = c(1.2,1,1,1,1),common.legend = T,legend = "bottom")
      
         }
    print(i)
   print(p_1)
    p[[j]]=p_1
    i=i+5
    j=j+1
    print(i)
  }
  p_all<-ggarrange(p[[1]],p[[2]],p[[3]],nrow=3,ncol=1,common.legend = T,legend="bottom",heights = c(1,1,1.8))
  return(p_all)  
}

plotSubSig_flip<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name,cols_x="black"){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group,col=I(cols_x)))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"]="NS"
  if(dim(stat.test)[1]>0){
    p1=p1 + coord_flip()+annotate("text",1:length(stat.test$p.adj.signif),max(stat.test$y.position),label=stat.test$p.adj.signif,col="red")
    
  }
  return(p1)
}


plotSubSig_pvalue<- function(x.data,xlab1,ylab1,sub1,title.x,cols1="",xlab.name,ylab.name){
  x.data<-x.data[,c(xlab1,ylab1,sub1)]
  colnames(x.data)<-c("x","y","group")
  p1=ggplot(data = x.data,aes(x=x,y=y))+geom_boxplot(aes(fill=group))+xlab(xlab.name)+ylab(ylab.name)+theme_hd()+
    ggtitle(title.x)+theme(axis.text.x = element_text(angle = 75,hjust = 1))
  if(cols1!=""){
    p1=p1+scale_fill_manual(values=cols1)
  }
  print(sub1)
  stat.test <- x.data %>%
    group_by( x) %>%
    wilcox_test(y ~group)  %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>% add_xy_position(x="x")
  stat.test=stat.test %>% filter(p.adj.signif!="ns")
  if(dim(stat.test)[1]>0){
    p1=p1 +  geom_signif(tip_length = 0.01,xmin=stat.test$xmin, xmax=stat.test$xmax, annotations=stat.test$p.adj, y_position=stat.test$y.position,col="black",textsize = 4.8,fontface="bold") 
    
  }
  return(p1)
}

plot3utrFreq<-function(genes_infor,temp_gene){
  
  
  temp_gene_trans<-aggregate(.~ensembl_gene_id+ensembl_transcript_id,data=temp_gene,max)
  
  temp_trans<-unique(unique(temp_gene[,1:2]))
  trans_count=table(temp_trans$ensembl_gene_id)
  
  temp_gene_uniq<-unique(temp_gene_trans[,-2])
  utr3_count<-table(temp_gene_uniq$ensembl_gene_id)
  genes_infor$utr3_count<-as.numeric(utr3_count[genes_infor$ensembl_gene_id])
  genes_infor$trans_count<-as.numeric(trans_count[genes_infor$ensembl_gene_id])
  genes_infor$utr3_freq<-genes_infor$utr3_count/genes_infor$trans_count
  
  
  temp_data=c()
  for(i in 0:10){
    temp_data<-rbind(temp_data,cbind(transript_count=paste(">",i),genes_infor[genes_infor$trans_count>i,]))
  }  
  temp_data[,"transript_count"]<-factor(temp_data[,"transript_count"],levels = paste(">",0:10))
  return(temp_data)
}



plot_histoligical<-function(genes_infor,col_x,comps,ylab){
  overlap.data=rbind(cbind(cancer="Carcinoma",genes_infor[which(genes_infor$Carcinoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Leukemia",genes_infor[which(genes_infor$Leukemia!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Lymphoma",genes_infor[which(genes_infor$Lymphoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Sarcoma",genes_infor[which(genes_infor$Sarcoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Myeloma",genes_infor[which(genes_infor$Myeloma!="" & genes_infor$Gene.Type=="TSG"),col_x]))
  overlap.data<-as.data.frame(overlap.data)
  overlap.data$cancer<-factor(overlap.data$cancer,levels=c("Myeloma","Carcinoma","Sarcoma","Leukemia","Lymphoma"))
  colnames(overlap.data)<-c("cancer_types",col_x)
  overlap.data[,col_x]<-as.numeric(overlap.data[,col_x])
  
  
  comp<-comps
    #list(c("Carcinoma","Leukemia"),c("Carcinoma","Lymphoma"),c("Carcinoma","Sarcoma"),c("Leukemia","Lymphoma"),c("Leukemia","Sarcoma"),c("Lymphoma","Sarcoma"),c("Myeloma","Leukemia"),c("Myeloma","Lymphoma"),c("Myeloma","Sarcoma"),c("Carcinoma","Myeloma"))
  #,test.args = c(alternative = "greater")
  px2_dif=ggplot(data = overlap.data,aes(x=cancer_types,y=get(col_x),fill=cancer_types))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold") +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"),axis.title.x =element_blank(),axis.text.x = element_text(angle =90))# 30,vjust=0.8
  px2_dif=px2_dif#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold",angle=90)
  print(px2_dif)
  genes_infor$cancer_count2=as.character(genes_infor$cancer_count2)
  genes_infor$cancer_count2[genes_infor$cancer_count2 %in% c("2",">2")]=">1"
  genes_infor$cancer_count2<-factor(genes_infor$cancer_count2,levels = c("1",">1"))
  px3_dif=ggplot(data = genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2 %in% c("1",">1"), ],aes(x=cancer_count2,y=get(col_x),fill=cancer_count2))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = list(c("1",">1")), test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold",test.args = c(alternative = "less")) +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"))
  px3_dif=px3_dif+xlab("Number of cancer type")+scale_fill_manual(values = c("#8EE5EE", "#20B2AA"))
  print(px3_dif)
  #+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  
  pall5<-list(p1=px2_dif,p2=px3_dif)
  return(pall5)
}



CutFile <- function(Input.name, Output.name, Annotation.name, Exception.symbols = "!Ser" ){
  Exception.symbols=toupper(Exception.symbols)
  if (file.exists(Output.name)){
    cat("Cut-File is already exited!\n")
    return(Output.name)
  }
  Input <- file(Input.name, "r")
  cat("",file=Output.name,append = FALSE)
  cat("",file=Annotation.name,append = FALSE)
  N.Char <- nchar(Exception.symbols)
  while (TRUE) {
    Line <- readLines(Input, n = 1)
    if (length(Line) == 0 ) {
      break
    }
    if ( toupper(substr(Line, start = 1, stop = N.Char)) == Exception.symbols ){
      cat(Line,"\n",file=Annotation.name,append = TRUE)
      next
    }
    if (!nzchar(Line)){
      next
    }
    cat(Line,"\n",file=Output.name,append = TRUE)
  }
  close(Input)
  return(Output.name)
}

getComp<-function(x){
  re_list<-list()
  for(i in 1:(length(x)-1)){
    for(j in i:(length(x)-1)){
      re_list[[paste(i,j)]]=c(x[i],x[j+1])
    }
    
  }
  return(re_list)
}


plot_histoligical<-function(genes_infor,col_x,comps,ylab){
  overlap.data=rbind(cbind(cancer="Carcinoma",genes_infor[which(genes_infor$Carcinoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Leukemia",genes_infor[which(genes_infor$Leukemia!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Lymphoma",genes_infor[which(genes_infor$Lymphoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Sarcoma",genes_infor[which(genes_infor$Sarcoma!="" & genes_infor$Gene.Type=="TSG"),col_x]),
                     cbind(cancer="Myeloma",genes_infor[which(genes_infor$Myeloma!="" & genes_infor$Gene.Type=="TSG"),col_x]))
  overlap.data<-as.data.frame(overlap.data)
  overlap.data$cancer<-factor(overlap.data$cancer,levels=c("Myeloma","Carcinoma","Sarcoma","Leukemia","Lymphoma"))
  colnames(overlap.data)<-c("cancer_types",col_x)
  overlap.data[,col_x]<-as.numeric(overlap.data[,col_x])
  
  
  comp<-comps
  #list(c("Carcinoma","Leukemia"),c("Carcinoma","Lymphoma"),c("Carcinoma","Sarcoma"),c("Leukemia","Lymphoma"),c("Leukemia","Sarcoma"),c("Lymphoma","Sarcoma"),c("Myeloma","Leukemia"),c("Myeloma","Lymphoma"),c("Myeloma","Sarcoma"),c("Carcinoma","Myeloma"))
  #,test.args = c(alternative = "greater")
  px2_dif=ggplot(data = overlap.data,aes(x=cancer_types,y=get(col_x),fill=cancer_types))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = comp, test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold") +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"),axis.title.x =element_blank(),axis.text.x = element_text(angle =90))# 30,vjust=0.8
  px2_dif=px2_dif#+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold",angle=90)
  print(px2_dif)
  genes_infor$cancer_count2=as.character(genes_infor$cancer_count2)
  genes_infor$cancer_count2[genes_infor$cancer_count2 %in% c("2",">2")]=">1"
  genes_infor$cancer_count2<-factor(genes_infor$cancer_count2,levels = c("1",">1"))
  px3_dif=ggplot(data = genes_infor[genes_infor$Gene.Type=="TSG" & genes_infor$cancer_count2 %in% c("1",">1"), ],aes(x=cancer_count2,y=get(col_x),fill=cancer_count2))+geom_boxplot(show.legend = FALSE,width=0.5,size=1)+
    geom_signif(comparisons = list(c("1",">1")), test='wilcox.test', map_signif_level=T,step_increase = 0.1,col="black",textsize = 4.8,fontface="bold",test.args = c(alternative = "less")) +ylab(ylab)+theme_hd()+
    ggtitle("TSGs")+ theme(plot.title = element_text(face = "bold"))
  px3_dif=px3_dif+xlab("Number of cancer type")+scale_fill_manual(values = c("#8EE5EE", "#20B2AA"))
  print(px3_dif)
  #+stat_summary(fun.data =  countFunction, geom="text", size = 4.8,col="black",fontface="bold")
  
  pall5<-list(p1=px2_dif,p2=px3_dif)
  return(pall5)
}

plot_classified_hist<-function(draws,cut_off,lables,x_lab){
  dens <- density(draws)
  dd <- with(dens,data.frame(x,y))
  cols_x=brewer.pal(n = 8, name = "Set2")
  # cols_x=lables
  s1=ggplot(data=dd,aes(x,y))+geom_line()+
    geom_ribbon(data=subset(dd,x<=cut_off[2]),aes(ymax=y,fill=I(cols_x[1])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[2] & x<=cut_off[3]),aes(ymax=y,fill=I(cols_x[2])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[3] & x<=cut_off[4]),aes(ymax=y,fill=I(cols_x[3])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd, x>cut_off[4] & x<=cut_off[5]),aes(ymax=y,fill=I(cols_x[4])),ymin=0,
                colour=NA)+
    geom_ribbon(data=subset(dd,x>cut_off[5]),aes(ymax=y,fill=I(cols_x[5])),ymin=0,
                colour=NA) + labs(fill = "Groups")+ylab("Density")+xlab(x_lab)
  
  print(s1)
}


plot_scatter<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,file_name){
  jpeg(file_name,width = 1300,height = 1200,res = 120,quality = 100)
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.1),color="cancer",size=0.6,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,label.y =label_y_pos)+facet_wrap(~cancer,nrow=3)  + theme(legend.position = "none",strip.text = element_text(size=25),axis.title =element_text(size=20))
  print(sp)
  dev.off()
}

plot_scatter_3<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", formatC(pvalue, format = "e", digits = 2),sep="")
  sp=sp+annotate(geom="text", x=label_x_pos-0.05, y=label_y_pos,size=4, label=data_text_1, color="black")
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}

plot_scatter_4<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", pvalue,sep="")
  sp=sp+annotate(geom="text", x=label_x_pos+0.1, y=label_y_pos,size=4, label=data_text_1, color="black")
  
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}


plot_scatter_5<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  all_data_tsg$Gene.Type=factor(all_data_tsg$Gene.Type,levels=c("Non-Cancer","TSG"))
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=1,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",label.x = label_x_pos,aes(color = get(col_x)),fontface="bold",output.type="text")  + theme(strip.text = element_text(size=25),axis.title =element_text(size=15))
  sc1<-cor.test(all_data_tsg[,y_lab],all_data_tsg[,x_lab])
  pvalue=sc1$p.value
  if(pvalue==0){pvalue="2.2e-16"}
  data_text_1= paste("All genes: R=",round(sc1$estimate,3),"   p < ", pvalue,sep="")
  sp=sp+facet_wrap(~cancer,nrow=5)+theme_hd()
  
  return(sp+scale_color_manual(values=c("#1874CD","#EE2C2C")))
}



plot_scatter_2<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.8),color=col_x,size=1.5,shape=16,
                  add = "reg.line",  add.params = list(color = "red"),  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}


plot_scatter_2kendall<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.6),color=col_x,size=1,shape=16,
                  add = "reg.line",  add.params = list(color = "blue"),  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "kendall",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}






plot_scatter_6<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.4),color=col_x,size=2,
                  add = "reg.line", add.params = list(color = "red"), # Add regressin line#size=1.5,shape=16,
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}


plot_scatter_7<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.1),color=col_x,size=0.7,shape=16,
                  add = "reg.line",  # Add regressin line#
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "pearson",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}



plot_scatter_2speaman<-function(all_data_tsg,x_lab,y_lab,x_title,y_title,label_x_pos,label_y_pos,col_x){
  sp <- ggscatter(all_data_tsg, 
                  x = x_lab, y = y_lab,alpha=I(0.8),color=col_x,size=1.5,shape=16,
                  add = "reg.line",  # Add regressin line
                  conf.int = TRUE # Add confidence interval
  )+xlab(x_title)+ylab(y_title)+theme_hd()
  # Add correlation coefficient
  sp=sp + stat_cor(method = "spearman",fontface="bold",output.type="text")  + theme(strip.text = element_text(size=12),axis.title =element_text(size=12))
  return(sp)
}

